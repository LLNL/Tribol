// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "QuadratureMatrixTransfer.hpp"

#include <unordered_map>

#include "axom/slic.hpp"

#include "redecomp/utils/MPIUtility.hpp"
#include "redecomp/utils/MPIArray.hpp"
#include "redecomp/RedecompMesh.hpp"

namespace redecomp
{

QuadratureMatrixTransfer::QuadratureMatrixTransfer(
  const mfem::QuadratureSpace& parent_test_space,
  const mfem::QuadratureSpace& parent_trial_space,
  const mfem::QuadratureSpace& redecomp_test_space,
  const mfem::QuadratureSpace& redecomp_trial_space
)
: parent_test_space_ { parent_test_space },
  parent_trial_space_ { parent_trial_space },
  redecomp_test_space_ { redecomp_test_space },
  redecomp_trial_space_ { redecomp_trial_space },
  redecomp_test_mesh_ { getRedecompMesh(redecomp_test_space_) },
  redecomp_trial_mesh_ { getRedecompMesh(redecomp_trial_space_) }
{
  SLIC_ERROR_ROOT_IF(static_cast<const mfem::Mesh*>(&redecomp_test_mesh_.getParent()) != parent_test_space_.GetMesh(),
    "The parent test quadrature space mesh must be linked to the test Redecomp mesh.");
  SLIC_ERROR_ROOT_IF(static_cast<const mfem::Mesh*>(&redecomp_trial_mesh_.getParent()) != parent_trial_space_.GetMesh(),
    "The parent trial quadrature space mesh must be linked to the trial Redecomp mesh.");
  SLIC_ERROR_ROOT_IF(&redecomp_test_mesh_.getMPIUtility().MPIComm() != &redecomp_trial_mesh_.getMPIUtility().MPIComm(),
    "MPI Communicator must match in test and trial spaces.");

  // one value per element for initial implementation
  SLIC_ERROR_ROOT_IF(parent_test_space_.GetOrder() != 1 || parent_trial_space_.GetOrder() != 1 ||
    redecomp_test_space_.GetOrder() != 1 || redecomp_trial_space_.GetOrder() != 1, 
    "Only one quadrature point per element is currently supported.");
}

std::tuple<
  std::unique_ptr<mfem::HypreParMatrix>, 
  axom::Array<int>, 
  axom::Array<int>
> QuadratureMatrixTransfer::TransferToParallel(
  const mfem::SparseMatrix& src
) const
{
  // check sizing on src
  SLIC_ERROR_ROOT_IF(src.Height() != redecomp_test_space_.GetSize(),
    "The number of rows in src must match the size of the test QuadratureSpace.");
  SLIC_ERROR_ROOT_IF(src.Width() != redecomp_trial_space_.GetSize(),
    "The number of columns in src must match the size of the trial QuadratureSpace.");
  
  // we want a csr matrix.  make sure src is finalized
  SLIC_ERROR_ROOT_IF(!src.Finalized(), "src must be finalized first.");

  auto my_rank = getMPIUtility().MyRank();
  auto n_ranks = getMPIUtility().NRanks();

  // package the entries to send to other ranks
  MPIArray<double> recv_matrix_data(&getMPIUtility());
  MPIArray<int> recv_parent_local_row(&getMPIUtility());
  MPIArray<int> recv_parent_local_col(&getMPIUtility());
  MPIArray<int> recv_parent_col_ranks(&getMPIUtility());
  {
    MPIArray<double> send_matrix_data(&getMPIUtility());
    MPIArray<int> send_parent_local_row_offsets(&getMPIUtility());
    MPIArray<int> send_parent_local_col_offsets(&getMPIUtility());
    MPIArray<int> send_parent_col_ranks(&getMPIUtility());
    // guess size and preallocate
    auto approx_size = 2 * src.NumNonZeroElems() / n_ranks;
    for (int r{0}; r < n_ranks; ++r)
    {
      send_matrix_data[r].reserve(approx_size);
      send_parent_local_row_offsets[r].reserve(approx_size);
      send_parent_local_col_offsets[r].reserve(approx_size);
      send_parent_col_ranks[r].reserve(approx_size);
    }

    auto src_I = src.GetI();
    auto src_J = src.GetJ();
    auto src_data = src.GetData();

    auto& test_elem_offsets = redecomp_test_mesh_.getRedecompToParentElemOffsets();
    auto& test_ghost_elems = redecomp_test_mesh_.getRedecompToParentGhostElems();
    auto& trial_elem_offsets = redecomp_trial_mesh_.getRedecompToParentElemOffsets();

    // test space = row space
    auto test_parent_rank = 0;
    auto test_ghost_ct = 0;
    for (int i{0}; i < src.Height(); ++i)
    {
      const auto test_elem_id = i; // assumes 1 point per element
      while (test_elem_id == test_elem_offsets[test_parent_rank + 1])
      {
        ++test_parent_rank;
        test_ghost_ct = 0;
      }
      // skip unowned elements
      if (test_ghost_ct < test_ghost_elems[test_parent_rank].size() 
        && test_elem_id == test_ghost_elems[test_parent_rank][test_ghost_ct])
      {
        ++test_ghost_ct;
        continue;
      }
      // trial space = col space
      auto trial_parent_rank = 0;
      for (int Ij{src_I[i]}; Ij < src_I[i+1]; ++Ij)
      {
        const int j{src_J[Ij]};
        const auto trial_elem_id = j; // assumes 1 point per element
        while (trial_elem_id >= trial_elem_offsets[trial_parent_rank + 1])
        {
          ++trial_parent_rank;
        }
        // if the entries in J aren't ordered, we might need to go backwards too
        while (trial_parent_rank > 0 && trial_elem_id < trial_elem_offsets[trial_parent_rank])
        {
          --trial_parent_rank;
        }
        // we have everything now.  package it up into send arrays
        send_matrix_data[test_parent_rank].push_back(src_data[Ij]);
        send_parent_local_row_offsets[test_parent_rank].push_back(test_elem_id - test_elem_offsets[test_parent_rank]);
        send_parent_local_col_offsets[test_parent_rank].push_back(trial_elem_id - trial_elem_offsets[trial_parent_rank]);
        send_parent_col_ranks[test_parent_rank].push_back(trial_parent_rank);
      }
    }
    // do MPI communication
    recv_matrix_data.SendRecvArrayEach(send_matrix_data);
    // NOTE: these hold offsets now; need to convert these to element numbers
    recv_parent_local_row.SendRecvArrayEach(send_parent_local_row_offsets);
    // NOTE: these are offsets, but if this rank doesn't own the element, it can't convert it to an element number.
    // need to communicate the offsets to the owned rank, get the local element number, then move it back
    recv_parent_local_col.SendRecvArrayEach(send_parent_local_col_offsets);
    recv_parent_col_ranks.SendRecvArrayEach(send_parent_col_ranks);
    // convert offsets to element numbers
    for (int r{0}; r < n_ranks; ++r)
    {
      for (auto& row_elem : recv_parent_local_row[r])
      {
        row_elem = redecomp_test_mesh_.getParentToRedecompElems().first[r][row_elem];
      }
      // do per-rank communication to transform local_col from offsets to element ids
      auto send_col_offsets_by_rank = MPIArray<int>(&getMPIUtility());
      // guess size and preallocate
      auto approx_rank_size = 2 * recv_parent_local_col[r].size() / n_ranks;
      for (auto& send_col_offsets : send_col_offsets_by_rank)
      {
        send_col_offsets.reserve(approx_rank_size);
      }
      for (int i{0}; i < recv_parent_local_col[r].size(); ++i)
      {
        send_col_offsets_by_rank[recv_parent_col_ranks[r][i]].push_back(recv_parent_local_col[r][i]);
      }
      auto recv_col_offsets_by_rank = MPIArray<int>(&getMPIUtility());
      recv_col_offsets_by_rank.SendRecvArrayEach(send_col_offsets_by_rank);
      // convert (now on-rank) offsets to element ids
      for (int trial_r{0}; trial_r < n_ranks; ++trial_r)
      {
        for (auto& recv_col_offset : recv_col_offsets_by_rank[trial_r])
        {
          recv_col_offset = redecomp_trial_mesh_.getParentToRedecompElems().first[trial_r][recv_col_offset];
        }
      }
      // send back to original rank
      send_col_offsets_by_rank.SendRecvArrayEach(recv_col_offsets_by_rank);
      // fill recv_parent_local_col with local element ids
      axom::Array<int> per_rank_ct(n_ranks, n_ranks);
      for (int i{0}; i < recv_parent_local_col[r].size(); ++i)
      {
        auto col_rank = recv_parent_col_ranks[r][i];
        recv_parent_local_col[r][i] = send_col_offsets_by_rank[col_rank][per_rank_ct[col_rank]++];
      }
    }
  }

  // Now we need to build the CSR data for the hypre diag and offd matrices.  diag holds data with both rows and cols
  // on-rank and offd holds data with rows on-rank and cols off-rank.
  auto row_starts = BuildParentElementRankOffsets(redecomp_test_mesh_);
  auto col_starts = BuildParentElementRankOffsets(redecomp_trial_mesh_);
  // size i_diag and i_offd, number of columms per row of matrix, and build a global dof to local offd column map
  auto diag_nnz = 0;
  auto offd_nnz = 0;
  // NOTE: these will be owned by the HypreParMatrix. using mfem::Memory, the data will not be deleted upon destruction.
  mfem::Memory<HYPRE_Int> i_diag(redecomp_test_space_.GetSize() + 1);
  mfem::Memory<HYPRE_Int> i_offd(redecomp_test_space_.GetSize() + 1);
  for (int i{0}; i < redecomp_test_space_.GetSize() + 1; ++i)
  {
    i_diag[i] = 0;
    i_offd[i] = 0;
  }
  std::unordered_map<HYPRE_BigInt, int> cmap_j_offd;
  for (int r{0}; r < n_ranks; ++r)
  {
    for (int i{0}; i < recv_parent_col_ranks[r].size(); ++i)
    {
      if (recv_parent_col_ranks[r][i] == my_rank)
      {
        ++diag_nnz;
        ++i_diag[recv_parent_local_row[r][i] + 1];
      }
      else
      {
        ++offd_nnz;
        ++i_offd[recv_parent_local_row[r][i] + 1];
        auto cmap_it = cmap_j_offd.insert(std::make_pair(
          col_starts[r] + recv_parent_local_col[r][i], cmap_j_offd.size()));
        recv_parent_local_col[r][i] = cmap_it.first->second;
      }
    }
  }
  // sum up i_diag and i_offd (completing their definition)
  for (int i{0}; i < redecomp_test_space_.GetSize(); ++i)
  {
    i_diag[i+1] += i_diag[i];
    i_offd[i+1] += i_offd[i];
  }
  SLIC_ASSERT(i_diag[redecomp_test_space_.GetSize()] == diag_nnz);
  SLIC_ASSERT(i_offd[redecomp_test_space_.GetSize()] == offd_nnz);
  // create cmap, offd column map from local to global element IDs
  mfem::Memory<HYPRE_BigInt> cmap(cmap_j_offd.size());
  for (auto& cmap_val : cmap_j_offd)
  {
    cmap[cmap_val.second] = cmap_val.first;
  }
  // compute j_diag and j_offd, the columns for each row (row offsets given by i).  will not be sorted.
  axom::Array<int> i_diag_ct(redecomp_test_space_.GetSize(), redecomp_test_space_.GetSize());
  axom::Array<int> i_offd_ct(redecomp_test_space_.GetSize(), redecomp_test_space_.GetSize());
  mfem::Memory<HYPRE_Int> j_diag(diag_nnz);
  mfem::Memory<double> data_diag(diag_nnz);
  mfem::Memory<HYPRE_Int> j_offd(offd_nnz);
  mfem::Memory<double> data_offd(offd_nnz);
  for (int r{0}; r < n_ranks; ++r)
  {
    for (int i{0}; i < recv_parent_col_ranks[r].size(); ++i)
    {
      auto curr_row = recv_parent_local_row[r][i];
      auto curr_col = recv_parent_local_col[r][i];
      if (recv_parent_col_ranks[r][i] == my_rank)
      {
        j_diag[i_diag[curr_row] + i_diag_ct[curr_row]] = curr_col;
        data_diag[i_diag[curr_row] + i_diag_ct[curr_row]] = recv_matrix_data[r][i];
        ++i_diag_ct[curr_row];
      }
      else
      {
        j_offd[i_offd[curr_row] + i_offd_ct[curr_row]] = curr_col;
        data_offd[i_offd[curr_row] + i_offd_ct[curr_row]] = recv_matrix_data[r][i];
        ++i_offd_ct[curr_row];
      }
    }
  }
  // resize row_starts and col_starts if needed
  if (HYPRE_AssumedPartitionCheck())
  {
    row_starts[0] = row_starts[my_rank];
    row_starts[1] = row_starts[my_rank + 1];
    if (row_starts.size() < 3)
    {
      auto old_size = row_starts.size();
      row_starts.resize(3);
      row_starts[2] = row_starts[old_size - 1];
    }
    else
    {
      row_starts[2] = row_starts.back();
      row_starts.resize(3);
      row_starts.shrink();
    }
    col_starts[0] = col_starts[my_rank];
    col_starts[1] = col_starts[my_rank + 1];
    if (col_starts.size() < 3)
    {
      auto old_size = col_starts.size();
      col_starts.resize(3);
      col_starts[2] = col_starts[old_size - 1];
    }
    else
    {
      col_starts[2] = col_starts.back();
      col_starts.resize(3);
      col_starts.shrink();
    }
  }
  // create hypre par matrix
  auto J_full = std::make_unique<mfem::HypreParMatrix>(getMPIUtility().MPIComm(),
    row_starts.back(), col_starts.back(), row_starts.data(), col_starts.data(),
    i_diag, j_diag, data_diag, i_offd, j_offd, data_offd, cmap_j_offd.size(), cmap
  );
  
  return {std::move(J_full), row_starts, col_starts};
}

const RedecompMesh& QuadratureMatrixTransfer::getRedecompMesh(const mfem::QuadratureSpace& quadrature_space)
{
  auto redecomp_mesh = dynamic_cast<const RedecompMesh*>(quadrature_space.GetMesh());
  SLIC_ERROR_ROOT_IF(redecomp_mesh == nullptr, "The quadrature space must have a Redecomp mesh.");
  return *redecomp_mesh;
}

axom::Array<HYPRE_BigInt> QuadratureMatrixTransfer::BuildParentElementRankOffsets(const RedecompMesh& redecomp_mesh)
{
  auto n_ranks = redecomp_mesh.getMPIUtility().NRanks();
  auto my_rank = redecomp_mesh.getMPIUtility().MyRank();
  axom::Array<HYPRE_BigInt> elem_offsets(n_ranks + 1, n_ranks + 1);
  elem_offsets[my_rank + 1] = redecomp_mesh.getParent().GetNE();
  redecomp_mesh.getMPIUtility().Allreduce(&elem_offsets, MPI_SUM);
  for (int r{2}; r < elem_offsets.size(); ++r)
  {
    elem_offsets[r] += elem_offsets[r - 1];
  }
  return elem_offsets;
}

const MPIUtility& QuadratureMatrixTransfer::getMPIUtility() const
{
  return redecomp_test_mesh_.getMPIUtility();
}

} // end namespace redecomp
