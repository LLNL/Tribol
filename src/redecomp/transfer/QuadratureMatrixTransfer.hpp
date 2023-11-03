// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_REDECOMP_QUADRATUREMATRIXTRANSFER_HPP_
#define SRC_REDECOMP_QUADRATUREMATRIXTRANSFER_HPP_

#include <memory>

#include "mfem.hpp"

#include "axom/core.hpp"

namespace redecomp
{

class MPIUtility;
class RedecompMesh;

/**
 * @brief Transfers sparse matrices on QuadratureSpaces from RedecompMesh to mfem::ParMesh
 *
 * This class enables transferring an mfem::SparseMatrix whose rows and columns are linked to mfem::QuadratureSpaces on
 * a RedecompMesh to a mfem::HypreParMatrix whose rows and columns are linked to mfem::QuadratureSpaces on an
 * mfem::ParMesh.
 */
class QuadratureMatrixTransfer
{
public:
  /**
   * @brief Construct a new QuadratureMatrixTransfer object
   * 
   * @param parent_test_space Test quadrature space on parent mfem::ParMesh
   * @param parent_trial_space Trial quadrature space on parent mfem::ParMesh
   * @param redecomp_test_space Test quadrature space on RedecompMesh
   * @param redecomp_trial_space Trial quadrature space on RedecompMesh
   */
  QuadratureMatrixTransfer(
    const mfem::QuadratureSpace& parent_test_space,
    const mfem::QuadratureSpace& parent_trial_space,
    const mfem::QuadratureSpace& redecomp_test_space,
    const mfem::QuadratureSpace& redecomp_trial_space
  );

  /**
   * @brief Transfers RedecompMesh sparse matrix to parent mfem::ParMesh hypre par matrix
   *
   * @param src Sparse matrix whose rows and columns are linked to test and trial QuadratureSpaces on the RedecompMesh
   * @return Tuple holding 1) owned pointer to mfem::HypreParMatrix on the parent mesh, 2) axom::Array of row starts,
   * and 3) axom::Array of column starts
   *
   * @note The HypreParMatrix does not own the row starts or column starts.  The user must ensure the lifetime of the
   * row starts and column starts matches the lifetime of the HypreParMatrix.
   */
  std::tuple<std::unique_ptr<mfem::HypreParMatrix>, axom::Array<int>, axom::Array<int>> TransferToParallel(
    const mfem::SparseMatrix& src
  ) const;

private:
  /**
   * @brief Get the RedecompMesh from the quadrature space
   *
   * This method verifies the quadrature space does have a RedecompMesh before returning. Quadrature spaces with no
   * RedecompMesh trigger a SLIC_ERROR macro.
   * 
   * @param quadrature_space Quadrature space linked to a RedecompMesh
   * @return const RedecompMesh& Redecomp mesh associated with quadrature space
   */
  static const RedecompMesh& getRedecompMesh(const mfem::QuadratureSpace& quadrature_space);

  /**
   * @brief Build vector of parent element offsets for each MPI rank
   * 
   * @param redecomp_mesh RedecompMesh whose parent linked mesh the element offset vector is desired
   * @return axom::Array<HYPRE_BigInt> Array of element offsets for each MPI rank
   */
  static axom::Array<HYPRE_BigInt> BuildParentElementRankOffsets(const RedecompMesh& redecomp_mesh);

  /**
   * @brief Returns MPIUtility pointer for the MatrixTransfer object
   * 
   * @return const MPIUtility* 
   */
  const MPIUtility& getMPIUtility() const;

  /**
   * @brief Test quadrature space on the parent mfem::ParMesh
   */
  const mfem::QuadratureSpace& parent_test_space_;

  /**
   * @brief Trial quadrature space on the parent mfem::ParMesh
   */
  const mfem::QuadratureSpace& parent_trial_space_;

  /**
   * @brief Test quadrature space on the RedecompMesh
   */
  const mfem::QuadratureSpace& redecomp_test_space_;

  /**
   * @brief Trial quadrature space on the RedecompMesh
   */
  const mfem::QuadratureSpace& redecomp_trial_space_;

  /**
   * @brief Test RedecompMesh
   */
  const RedecompMesh& redecomp_test_mesh_;

  /**
   * @brief Trial RedecompMesh
   */
  const RedecompMesh& redecomp_trial_mesh_;

};

} // end namespace redecomp

#endif /* SRC_REDECOMP_QUADRATUREMATRIXTRANSFER_HPP_ */
