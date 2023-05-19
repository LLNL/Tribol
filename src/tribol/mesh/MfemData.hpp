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

// dual fields only go from redecomp to submesh
class SubmeshRedecompTransfer
{
public:
  SubmeshRedecompTransfer(
    redecomp::RedecompMesh& redecomp,
    const mfem::ParFiniteElementSpace& submesh_fes
  );
  mfem::GridFunction SubmeshToRedecomp(const mfem::ParGridFunction& src) const;
  mfem::ParGridFunction RedecompToSubmesh(const mfem::GridFunction& src) const;
  void RedecompToSubmesh(
    const mfem::GridFunction& src,
    mfem::ParGridFunction& dst
  ) const;
  mfem::ParGridFunction EmptySubmeshGridFn() const
  {
    return mfem::ParGridFunction(&submesh_fes_);
  }
  const mfem::ParSubMesh& GetSubmesh() const 
  { 
    return static_cast<const mfem::ParSubMesh&>(*submesh_fes_.GetParMesh());
  }
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);
private:
  mutable mfem::ParFiniteElementSpace submesh_fes_;
  mutable std::unique_ptr<mfem::FiniteElementSpace> redecomp_fes_;
  const redecomp::RedecompTransfer redecomp_xfer_;
};

class ParentRedecompTransfer
{
public:
  ParentRedecompTransfer(
    redecomp::RedecompMesh& redecomp,
    mfem::ParSubMesh& submesh,
    const mfem::ParFiniteElementSpace& parent_fes
  );
  mfem::GridFunction ParentToRedecomp(const mfem::ParGridFunction& src) const;
  mfem::ParGridFunction RedecompToParent(const mfem::GridFunction& src) const;
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);
private:
  mutable mfem::ParFiniteElementSpace parent_fes_;
  SubmeshRedecompTransfer redecomp_xfer_;
  mutable mfem::ParGridFunction submesh_gridfn_;
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
  std::vector<const real*> GetFieldPtrs() const;
  static std::vector<real*> GetFieldPtrs(mfem::GridFunction& redecomp_fn);
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
  std::vector<const real*> GetFieldPtrs() const;
  static std::vector<real*> GetFieldPtrs(mfem::GridFunction& redecomp_fn);
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
    const std::set<integer>& attributes_2,
    std::unique_ptr<mfem::FiniteElementCollection> dual_fec = nullptr,
    integer dual_vdim = 0
  );
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
  std::vector<const real*> GetCoordsPtrs() const 
  { 
    return coords_.GetFieldPtrs();
  }
  std::vector<real*> GetResponsePtrs()
  {
    return PrimalField::GetFieldPtrs(response_gridfn_);
  }
  mfem::ParGridFunction GetParentResponse() const
  {
    return GetPrimalTransfer().RedecompToParent(response_gridfn_);
  }
  const ParentRedecompTransfer& GetPrimalTransfer() const
  { 
    return GetUpdateData().primal_xfer_;
  }
  void SetParentVelocity(const mfem::ParGridFunction& velocity);
  bool HasVelocity() const { return velocity_ != nullptr; }
  std::vector<const real*> GetVelocityPtrs() const
  {
    return velocity_->GetFieldPtrs();
  }
  std::vector<const real*> GetPressurePtrs() const
  {
    return pressure_->GetFieldPtrs();
  }
  std::vector<real*> GetGapPtrs()
  {
    return PressureField::GetFieldPtrs(*gap_gridfn_);
  }
private:
  struct UpdateData
  {
    UpdateData(
      mfem::ParSubMesh& submesh,
      const std::set<integer>& attributes_1,
      const std::set<integer>& attributes_2,
      integer num_verts_per_elem,
      const mfem::ParFiniteElementSpace& parent_fes,
      const mfem::ParFiniteElementSpace* redecomp_fes
    );
    redecomp::RedecompMesh redecomp_;
    ParentRedecompTransfer primal_xfer_;
    std::unique_ptr<SubmeshRedecompTransfer> dual_xfer_;
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
  integer mesh_id_1_;
  integer mesh_id_2_;
  mfem::ParMesh& mesh_;
  std::set<integer> attributes_1_;
  std::set<integer> attributes_2_;
  mfem::ParSubMesh submesh_;
  PrimalField coords_;
  std::unique_ptr<mfem::ParGridFunction> pressure_gridfn_;
  std::unique_ptr<PressureField> pressure_;
  std::unique_ptr<mfem::GridFunction> gap_gridfn_;
  std::unique_ptr<const mfem::ParFiniteElementSpace> dual_fes_;
  integer dual_vdim_;
  mfem::GridFunction response_gridfn_;
  InterfaceElementType elem_type_;
  integer num_verts_per_elem_;
  std::unique_ptr<UpdateData> update_data_;
  std::unique_ptr<PrimalField> velocity_;

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

} // end namespace tribol

#endif /* SRC_MESH_MFEMDATA_HPP_ */
