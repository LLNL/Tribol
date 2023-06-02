// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MFEMDATA_HPP_
#define SRC_MESH_MFEMDATA_HPP_

#include <set>
#include <vector>

#include "mfem.hpp"

#include "axom/core.hpp"
#include "redecomp/redecomp.hpp"

#include "tribol/types.hpp"
#include "tribol/common/Parameters.hpp"

namespace tribol
{

class SubmeshLORTransfer
{
public:
  SubmeshLORTransfer(
    mfem::ParFiniteElementSpace& ho_submesh_fes,
    mfem::ParMesh& lor_submesh
  );
  void TransferToLORGridFn(
    const mfem::ParGridFunction& ho_src
  );
  void TransferFromLORGridFn(
    mfem::ParGridFunction& ho_dst
  ) const;
  mfem::ParGridFunction& GetLORGridFn() { return lor_gridfn_; }
  const mfem::ParGridFunction& GetLORGridFn() const { return lor_gridfn_; }
private:
  static mfem::ParGridFunction CreateLORGridFunction(
    mfem::ParMesh& lor_mesh,
    std::unique_ptr<mfem::FiniteElementCollection> lor_fec,
    integer vdim
  );
  mfem::ParGridFunction lor_gridfn_;
  mutable mfem::L2ProjectionGridTransfer lor_xfer_;
};

// dual fields only go from redecomp to submesh
class SubmeshRedecompTransfer
{
public:
  SubmeshRedecompTransfer(
    mfem::ParFiniteElementSpace& submesh_fes,
    SubmeshLORTransfer* lor_xfer,
    redecomp::RedecompMesh& redecomp
  );
  mfem::GridFunction SubmeshToRedecomp(const mfem::ParGridFunction& src) const;
  mfem::ParGridFunction RedecompToSubmesh(const mfem::GridFunction& src) const;
  void RedecompToSubmesh(
    const mfem::GridFunction& src,
    mfem::ParGridFunction& dst
  ) const;
  const mfem::ParSubMesh& GetSubmesh() const 
  { 
    return static_cast<const mfem::ParSubMesh&>(*submesh_fes_.GetParMesh());
  }
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);
private:
  static std::unique_ptr<mfem::FiniteElementSpace> CreateRedecompFESpace(
    redecomp::RedecompMesh& redecomp,
    mfem::ParFiniteElementSpace& submesh_fes
  );
  mfem::ParFiniteElementSpace& submesh_fes_;
  mutable std::unique_ptr<mfem::FiniteElementSpace> redecomp_fes_;
  SubmeshLORTransfer* lor_xfer_;
  const redecomp::RedecompTransfer redecomp_xfer_;
};

class ParentRedecompTransfer
{
public:
  ParentRedecompTransfer(
    const mfem::ParFiniteElementSpace& parent_fes,
    mfem::ParGridFunction& submesh_gridfn,
    SubmeshLORTransfer* lor_xfer,
    redecomp::RedecompMesh& redecomp
  );
  mfem::GridFunction ParentToRedecomp(const mfem::ParGridFunction& src) const;
  mfem::ParGridFunction RedecompToParent(const mfem::GridFunction& src) const;
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_gridfn_.ParFESpace();
  }
private:
  mutable mfem::ParFiniteElementSpace parent_fes_;
  mfem::ParGridFunction& submesh_gridfn_;
  SubmeshRedecompTransfer redecomp_xfer_;
};

class PrimalField
{
public:
  PrimalField(const mfem::ParGridFunction& parent);
  void SetParentField(const mfem::ParGridFunction& parent);
  void UpdateField(const ParentRedecompTransfer& xfer);
  const mfem::ParGridFunction& GetParentField() const
  {
    return parent_;
  }
  mfem::GridFunction& GetRedecompField() { return GetUpdateData().redecomp_; }
  const mfem::GridFunction& GetRedecompField() const
  {
    return GetUpdateData().redecomp_;
  }
  std::vector<const real*> GetRedecompFieldPtrs() const;
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn);
private:
  struct UpdateData
  {
    UpdateData(
      const ParentRedecompTransfer& xfer,
      const mfem::ParGridFunction& parent
    );
    const ParentRedecompTransfer& xfer_;
    mfem::GridFunction redecomp_;
  };
  UpdateData& GetUpdateData();
  const UpdateData& GetUpdateData() const;
  std::reference_wrapper<const mfem::ParGridFunction> parent_;
  std::unique_ptr<UpdateData> update_data_;
};

class PressureField
{
public:
  PressureField(const mfem::ParGridFunction& submesh);
  void SetSubmeshField(const mfem::ParGridFunction& submesh);
  void UpdateField(const SubmeshRedecompTransfer& xfer);
  const mfem::ParGridFunction& GetSubmeshField() const
  {
    return submesh_;
  }
  mfem::GridFunction& GetRedecompField() { return GetUpdateData().redecomp_; }
  const mfem::GridFunction& GetRedecompField() const
  {
    return GetUpdateData().redecomp_;
  }
  std::vector<const real*> GetRedecompFieldPtrs() const;
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn);
private:
  struct UpdateData
  {
    UpdateData(
      const SubmeshRedecompTransfer& xfer,
      const mfem::ParGridFunction& submesh
    );
    const SubmeshRedecompTransfer& xfer_;
    mfem::GridFunction redecomp_;
  };
  UpdateData& GetUpdateData();
  const UpdateData& GetUpdateData() const;
  std::reference_wrapper<const mfem::ParGridFunction> submesh_;
  std::unique_ptr<UpdateData> update_data_;
};

class MfemMeshData
{
public:
  MfemMeshData(
    integer mesh_id_1,
    integer mesh_id_2,
    mfem::ParMesh& mesh,
    const mfem::ParGridFunction& current_coords,
    const std::set<integer>& attributes_1,
    const std::set<integer>& attributes_2
  );
  const mfem::ParGridFunction& GetParentCoords() const
  {
    return coords_.GetParentField();
  }
  void SetParentCoords(const mfem::ParGridFunction& current_coords);
  void UpdateMeshData();
  integer GetMesh1ID() const { return mesh_id_1_; }
  integer GetMesh2ID() const { return mesh_id_2_; }
  integer GetMesh1NE() const
  { 
    return GetUpdateData().conn_1_.size() / num_verts_per_elem_;
  }
  integer GetMesh2NE() const
  {
    return GetUpdateData().conn_2_.size() / num_verts_per_elem_;
  }
  integer GetNV() const { return GetUpdateData().redecomp_.GetNV(); }
  const IndexType* GetMesh1Conn() const
  {
    return GetUpdateData().conn_1_.data();
  }
  const IndexType* GetMesh2Conn() const
  {
    return GetUpdateData().conn_2_.data();
  }
  integer GetElemType() const { return elem_type_; }
  std::vector<const real*> GetRedecompCoordsPtrs() const 
  { 
    return coords_.GetRedecompFieldPtrs();
  }
  std::vector<real*> GetRedecompResponsePtrs()
  {
    return PrimalField::GetRedecompFieldPtrs(redecomp_response_);
  }
  const mfem::GridFunction& GetRedecompResponse() const
  {
    return redecomp_response_;
  }
  mfem::ParGridFunction GetParentResponse() const
  {
    return GetPrimalTransfer().RedecompToParent(redecomp_response_);
  }
  const ParentRedecompTransfer& GetPrimalTransfer() const
  { 
    return GetUpdateData().primal_xfer_;
  }
  void SetParentVelocity(const mfem::ParGridFunction& velocity);
  bool HasVelocity() const { return velocity_ != nullptr; }
  std::vector<const real*> GetRedecompVelocityPtrs() const
  {
    return velocity_->GetRedecompFieldPtrs();
  }
  const std::vector<integer>& GetElemMap1() const
  {
    return GetUpdateData().elem_map_1_;
  }
  const std::vector<integer>& GetElemMap2() const
  {
    return GetUpdateData().elem_map_2_;
  }
  mfem::ParSubMesh& GetSubmesh() { return submesh_; }
  mfem::ParMesh* GetLORSubmesh() { return lor_submesh_.get(); }
  redecomp::RedecompMesh& GetRedecompMesh()
  {
    return GetUpdateData().redecomp_;
  }
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return GetUpdateData().primal_xfer_.GetSubmeshFESpace();
  }
  void SetLowOrderRefinedFactor(integer lor_factor);
private:
  struct UpdateData
  {
    UpdateData(
      mfem::ParSubMesh& submesh,
      mfem::ParMesh* lor_submesh,
      const mfem::ParFiniteElementSpace& parent_fes,
      mfem::ParGridFunction& submesh_gridfn,
      SubmeshLORTransfer* lor_xfer,
      const std::set<integer>& attributes_1,
      const std::set<integer>& attributes_2,
      integer num_verts_per_elem
    );
    redecomp::RedecompMesh redecomp_;
    ParentRedecompTransfer primal_xfer_;
    axom::Array<integer, 2> conn_1_;
    axom::Array<integer, 2> conn_2_;
    std::vector<integer> elem_map_1_;
    std::vector<integer> elem_map_2_;
  private:
    void UpdateConnectivity(
      const std::set<integer>& attributes_1,
      const std::set<integer>& attributes_2,
      integer num_verts_per_elem
    );
  };
  UpdateData& GetUpdateData();
  const UpdateData& GetUpdateData() const;
  static mfem::ParSubMesh CreateSubmesh(
    const mfem::ParMesh& parent,
    const std::set<integer>& attributes_1,
    const std::set<integer>& attributes_2
  );
  static mfem::ParGridFunction CreateSubmeshGridFn(
    mfem::ParSubMesh& submesh,
    mfem::ParFiniteElementSpace& parent_fes
  );
  integer mesh_id_1_;
  integer mesh_id_2_;
  mfem::ParMesh& mesh_;
  std::set<integer> attributes_1_;
  std::set<integer> attributes_2_;
  mfem::ParSubMesh submesh_;
  PrimalField coords_;
  mfem::ParGridFunction submesh_xfer_gridfn_;
  integer lor_factor_;
  std::unique_ptr<mfem::ParMesh> lor_submesh_;
  std::unique_ptr<SubmeshLORTransfer> lor_xfer_;
  InterfaceElementType elem_type_;
  integer num_verts_per_elem_;
  std::unique_ptr<PrimalField> velocity_;
  std::unique_ptr<UpdateData> update_data_;
  mfem::GridFunction redecomp_response_;

  template <typename T>
  static T mergeContainers(T container_1, T container_2)
  {
    auto merged = container_1;
    merged.insert(container_2.begin(), container_2.end());
    return merged;
  }

  /**
   * @brief Converts a std::set to an mfem::Array
   * 
   * @tparam T type held in the set and array
   * @param orig original set
   * @return mfem::Array<T> output array holding entries in orig
   */
  template <typename T>
  static mfem::Array<T> arrayFromSet(std::set<T> orig)
  {
    auto array = mfem::Array<T>();
    array.Reserve(static_cast<int>(orig.size()));
    for (const auto& val : orig)
    {
      array.Append(val);
    }
    return array;
  }
};

class MfemDualData
{
public:
  MfemDualData(
    mfem::ParSubMesh& submesh,
    mfem::ParMesh* lor_submesh,
    std::unique_ptr<mfem::FiniteElementCollection> dual_fec,
    integer dual_vdim
  );
  void UpdateDualData(redecomp::RedecompMesh& redecomp);
  std::vector<const real*> GetRedecompPressurePtrs() const
  {
    return pressure_.GetRedecompFieldPtrs();
  }
  mfem::ParGridFunction& GetSubmeshPressure()
  {
    return submesh_pressure_;
  }
  std::vector<real*> GetRedecompGapPtrs()
  {
    return PressureField::GetRedecompFieldPtrs(redecomp_gap_);
  }
  const mfem::GridFunction& GetRedecompGap() const
  {
    return redecomp_gap_;
  }
  mfem::ParGridFunction GetSubmeshGap() const
  {
    return GetDualTransfer().RedecompToSubmesh(redecomp_gap_);
  }
  const SubmeshRedecompTransfer& GetDualTransfer() const
  {
    return GetUpdateData().dual_xfer_;
  }
private:
  struct UpdateData
  {
    UpdateData(
      mfem::ParFiniteElementSpace& submesh_fes,
      SubmeshLORTransfer* lor_xfer,
      redecomp::RedecompMesh& redecomp
    );
    SubmeshRedecompTransfer dual_xfer_;
  };
  UpdateData& GetUpdateData();
  const UpdateData& GetUpdateData() const;
  mfem::ParGridFunction submesh_pressure_;
  PressureField pressure_;
  std::unique_ptr<SubmeshLORTransfer> lor_xfer_;
  mfem::GridFunction redecomp_gap_;
  std::unique_ptr<UpdateData> update_data_;
};

} // end namespace tribol

#endif /* SRC_MESH_MFEMDATA_HPP_ */
