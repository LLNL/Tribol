// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MFEMDATA_HPP_
#define SRC_MESH_MFEMDATA_HPP_

#include "tribol/types.hpp"

#ifdef BUILD_REDECOMP

#include <set>
#include <vector>

#include "mfem.hpp"

#include "axom/core.hpp"
#include "redecomp/redecomp.hpp"

#include "tribol/common/Parameters.hpp"
#include "tribol/mesh/MethodCouplingData.hpp"

namespace tribol
{

/**
 * @brief Facilitates transfer of fields to/from parent-linked boundary submesh
 * (a higher-order mesh) to LOR mesh
 *
 * This class simplifies transferring 1) a primal field such as displacement and
 * velocity from a higher-order grid function on a parent-linked boundary
 * submesh to a low-order grid function on a LOR boundary mesh and 2) the
 * energetic conjugate to a primal field (such as nodal force, which is
 * conjugate to nodal displacement) from a low-order grid function on a LOR mesh
 * to a higher-order grid function on a parent-linked boundary submesh.
 *
 * Field data on the LOR mesh are stored internally in this class, and accessed
 * through the GetLORGridFn() method.
 */
class SubmeshLORTransfer
{
public:
  /**
   * @brief Construct a new SubmeshLORTransfer object
   *
   * @param submesh_fes Higher order finite element space on the parent-linked
   * boundary submesh
   * @param lor_mesh LOR mesh
   */
  SubmeshLORTransfer(
    mfem::ParFiniteElementSpace& submesh_fes,
    mfem::ParMesh& lor_mesh
  );

  /**
   * @brief Transfers data from a higher-order grid function on a parent-linked
   * submesh
   *
   * Data is transferred to the low-order grid function on the LOR mesh that can
   * be accessed using GetLORGridFn().
   *
   * @param submesh_src Source higher-order grid function on the parent-linked
   * boundary submesh
   */
  void TransferToLORGridFn(
    const mfem::ParGridFunction& submesh_src
  );

  /**
   * @brief Transfers data to a higher-order grid function on a parent-linked
   * boundary submesh
   *
   * Data must be stored in the low-order grid function on the LOR mesh accessed
   * using GetLORGridFn().
   *
   * @param submesh_dst Destination higher-order grid function on the
   * parent-linked boundary submesh
   */
  void TransferFromLORGridFn(
    mfem::ParGridFunction& submesh_dst
  ) const;

  /**
   * @brief Access the local low-order grid function on the LOR mesh
   * 
   * @return mfem::ParGridFunction& 
   */
  mfem::ParGridFunction& GetLORGridFn() { return lor_gridfn_; }

  /**
   * @brief Access the local low-order grid function on the LOR mesh
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetLORGridFn() const { return lor_gridfn_; }

private:
  /**
  * @brief Create low-order grid function on the LOR mesh
  * 
  * @param lor_mesh LOR mesh
  * @param lor_fec Finite element collection to apply to grid function
  * @param vdim Vector dimension of the grid function
  * @return mfem::ParGridFunction on lor_mesh, with lor_fec and vdim specified
  */
  static mfem::ParGridFunction CreateLORGridFunction(
    mfem::ParMesh& lor_mesh,
    std::unique_ptr<mfem::FiniteElementCollection> lor_fec,
    integer vdim
  );

  /**
   * @brief Local low-order grid function on the LOR mesh 
   */
  mfem::ParGridFunction lor_gridfn_;

  /**
   * @brief Low-order refined <-> higher-order coarse transfer object
   */
  mutable mfem::L2ProjectionGridTransfer lor_xfer_;
};

/**
 * @brief Facilitates transferring variables from the submesh to redecomp mesh
 * levels
 *
 * This class simplifies transferring field variables from/to a grid function on
 * an mfem::ParSubMesh to/from a grid function on a redecomp::RedecompMesh.  If
 * transferring to a low-order grid function on a LOR mesh is also required,
 * this class will perform that transfer as well.
 *
 * The hierarchy of transfer is: submesh <--> LOR mesh (optional) <--> redecomp
 * mesh
 *
 * @note This is used to transfer variables defined at the mfem::ParSubMesh
 * level (e.g. pressure and gap).
 */
class SubmeshRedecompTransfer
{
public:
  /**
   * @brief Construct a new SubmeshRedecompTransfer object
   *
   * @param submesh_fes Finite element space on the parent-linked boundary
   * submesh
   * @param submesh_lor_xfer Submesh to LOR grid function transfer object (if
   * using LOR; nullptr otherwise)
   * @param redecomp_mesh RedecompMesh of the redecomposed contact surface mesh
   */
  SubmeshRedecompTransfer(
    mfem::ParFiniteElementSpace& submesh_fes,
    SubmeshLORTransfer* submesh_lor_xfer,
    redecomp::RedecompMesh& redecomp_mesh
  );

  /**
   * @brief Transfer grid function on parent-linked boundary submesh to grid
   * function on redecomp mesh
   *
   * @param [in] submesh_src Grid function on parent-linked boundary submesh
   * @param [out] redecomp_dst Zero-valued grid function on redecomp mesh
   */
  void SubmeshToRedecomp(
    const mfem::ParGridFunction& submesh_src,
    mfem::GridFunction& redecomp_dst
  ) const;

  /**
   * @brief Transfer grid function on redecomp mesh to grid function on
   * parent-linked boundary submesh
   *
   * @param redecomp_src Grid function on redecomp mesh
   * @param submesh_dst Zero-valued grid function on parent-linked boundary
   * submesh
   */
  void RedecompToSubmesh(
    const mfem::GridFunction& redecomp_src,
    mfem::ParGridFunction& submesh_dst
  ) const;

  /**
   * @brief Get the parent-linked boundary submesh associated with the
   * SubmeshRedecompTransfer object
   *
   * @return const mfem::ParSubMesh& 
   */
  const mfem::ParSubMesh& GetSubmesh() const 
  { 
    return static_cast<const mfem::ParSubMesh&>(*submesh_fes_.GetParMesh());
  }

  /**
   * @brief Returns finite element space on the redecomp mesh associated with
   * this transfer object
   *
   * @return mfem::FiniteElementSpace& 
   */
  mfem::FiniteElementSpace& GetRedecompFESpace()
  {
    return *redecomp_fes_;
  }

private:
  /**
   * @brief Create a finite element space on the redecomp mesh
   *
   * @param redecomp_mesh RedecompMesh of the redecomposed contact surface mesh
   * @param submesh_fes Finite element space on the parent-linked boundary
   * submesh
   * @return std::unique_ptr<mfem::FiniteElementSpace> 
   */
  static std::unique_ptr<mfem::FiniteElementSpace> CreateRedecompFESpace(
    redecomp::RedecompMesh& redecomp_mesh,
    mfem::ParFiniteElementSpace& submesh_fes
  );

  /**
   * @brief Finite element space on the parent-linked boundary submesh 
   */
  mfem::ParFiniteElementSpace& submesh_fes_;

  /**
   * @brief Finite element space on the redecomp mesh
   */
  mutable std::unique_ptr<mfem::FiniteElementSpace> redecomp_fes_;

  /**
   * @brief Transfer object between low-order and higher-order grid functions
   */
  SubmeshLORTransfer* submesh_lor_xfer_;

  /**
   * @brief Transfer object between parent-linked boundary submesh and redecomp
   * mesh
   */
  const redecomp::RedecompTransfer redecomp_xfer_;
};

/**
 * @brief Facilitates transferring variables from the parent mesh to the
 * redecomp mesh levels
 *
 * This class simplifies transferring field variables from/to a grid function on
 * a parent mfem::ParMesh to/from a grid function on a redecomp::RedecompMesh.
 * If transferring to a low-order grid function on a LOR mesh is also required,
 * this class will perform that transfer as well.
 *
 * The hierarchy of transfer is:
 * parent mesh <--> submesh <--> LOR mesh (optional) <--> redecomp mesh
 *                 \---------------------------------------------------/
 *                handled through SubmeshRedecompTransfer member variable
 */
class ParentRedecompTransfer
{
public:
  /**
   * @brief Construct a new ParentRedecompTransfer object
   *
   * @param parent_fes Finite element space on the parent mesh
   * @param submesh_gridfn Grid function on the parent-linked boundary submesh
   * used to temporarily store variables being transferred
   * @param submesh_lor_xfer Submesh to LOR grid function transfer object (if
   * using LOR; nullptr otherwise)
   * @param redecomp_mesh RedecompMesh of the redecomposed contact surface mesh
   */
  ParentRedecompTransfer(
    const mfem::ParFiniteElementSpace& parent_fes,
    mfem::ParGridFunction& submesh_gridfn,
    SubmeshLORTransfer* submesh_lor_xfer,
    redecomp::RedecompMesh& redecomp_mesh
  );

  /**
   * @brief Transfer grid function on parent mesh to grid function on redecomp
   * mesh
   *
   * @param [in] parent_src Grid function on parent mesh
   * @param [out] redecomp_dst Zero-valued grid function on redecomp mesh
   */
  void ParentToRedecomp(
    const mfem::ParGridFunction& parent_src,
    mfem::GridFunction& redecomp_dst
  ) const;
  
  /**
   * @brief Transfer grid function on redecomp mesh to grid function on parent
   * mesh
   *
   * @param [in] redecomp_src Grid function on RedecompMesh
   * @param [out] parent_dst Zero-valued grid function on parent mesh
   */
  void RedecompToParent(
    const mfem::GridFunction& redecomp_src, 
    mfem::ParGridFunction& parent_dst
  ) const;

  /**
   * @brief Get the parent-linked boundary submesh finite element space
   * associated with this transfer object
   *
   * @return const mfem::ParFiniteElementSpace& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_gridfn_.ParFESpace();
  }

  /**
   * @brief Returns finite element space on the redecomp mesh associated with
   * this transfer object
   *
   * @return mfem::FiniteElementSpace& 
   */
  mfem::FiniteElementSpace& GetRedecompFESpace()
  {
    return submesh_redecomp_xfer_.GetRedecompFESpace();
  }

private:
  /**
   * @brief Finite element space on the parent mesh
   */
  mutable mfem::ParFiniteElementSpace parent_fes_;

  /**
   * @brief Grid function on the parent-linked boundary submesh
   */
  mfem::ParGridFunction& submesh_gridfn_;

  /**
   * @brief Object to transfer variables to/from the parent-linked boundary
   * submesh level from/to the redecomp level
   */
  SubmeshRedecompTransfer submesh_redecomp_xfer_;
};

/**
 * @brief Vector field variable that lives on the parent mesh
 *
 * This class stores a vector field variable defined on the parent mesh and
 * handles transferring field data to/from different mesh representations used
 * by Tribol.
 *
 * @note Example parent vector fields include displacement and velocity.
 */
class ParentField
{
public:
  /**
   * @brief Construct a new ParentField object
   * 
   * @param parent Grid function on the parent mesh
   */
  ParentField(const mfem::ParGridFunction& parent_gridfn);
  
  /**
   * @brief Set a new grid function on the parent mesh
   * 
   * @param parent Grid function on the parent mesh
   */
  void SetParentGridFn(const mfem::ParGridFunction& parent_gridfn);

  /**
   * @brief Set a new transfer object when the redecomp mesh has been updated
   *
   * @param xfer Updated parent mesh to redecomp mesh transfer object
   */
  void UpdateField(ParentRedecompTransfer& parent_redecomp_xfer);

  /**
   * @brief Get the parent grid function
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetParentGridFn() const
  {
    return parent_gridfn_;
  }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return mfem::GridFunction& 
   */
  mfem::GridFunction& GetRedecompGridFn() { return GetUpdateData().redecomp_gridfn_; }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompGridFn() const
  {
    return GetUpdateData().redecomp_gridfn_;
  }

  /**
   * @brief Get pointers to component arrays of the redecomp mesh grid function
   *
   * @return std::vector<const real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  std::vector<const real*> GetRedecompFieldPtrs() const;
  
  /**
   * @brief Get pointers to component arrays of the redecomp mesh grid function
   * 
   * @param redecomp_gridfn Redecomp mesh grid function
   * @return std::vector<real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_gridfn);

private:

  /**
   * @brief Creates and stores data that changes when the redecomp mesh is
   * updated
   */
  struct UpdateData
  {
    /**
     * @brief Construct a new UpdateData object
     * 
     * @param parent_redecomp_xfer Parent to redecomp field transfer object
     * @param parent_gridfn Grid function on the original, parent mesh
     */
    UpdateData(
      ParentRedecompTransfer& parent_redecomp_xfer,
      const mfem::ParGridFunction& parent_gridfn
    );

    /**
     * @brief Parent to redecomp field transfer object
     */
    const ParentRedecompTransfer& parent_redecomp_xfer_;

    /**
     * @brief Grid function values on the redecomp mesh
     */
    mfem::GridFunction redecomp_gridfn_;
  };

  /**
   * @brief Get the UpdateData object
   * 
   * @return UpdateData& 
   */
  UpdateData& GetUpdateData();

  /**
   * @brief Get the UpdateData object
   * 
   * @return const UpdateData& 
   */
  const UpdateData& GetUpdateData() const;

  /**
   * @brief Grid function on the parent mesh
   *
   * @note Stored as a reference wrapper so the reference can be updated
   */
  std::reference_wrapper<const mfem::ParGridFunction> parent_gridfn_;

  /**
   * @brief UpdateData object created upon call to UpdateField()
   */
  std::unique_ptr<UpdateData> update_data_;
};

/**
 * @brief Stores a pressure field variable that lives on the parent-linked
 * boundary submesh
 *
 * This class handles transferring pressure field data to/from representations
 * used by MFEM from/to representations used by Tribol.
 */
class PressureField
{
public:
  /**
   * @brief Construct a new PressureField object
   *
   * @param submesh_gridfn Grid function on the parent-linked boundary submesh
   */
  PressureField(const mfem::ParGridFunction& submesh_gridfn);

  /**
   * @brief Sets a new grid function on the parent-linked boundary submesh
   * 
   * @param submesh_gridfn Grid function on the parent-linked boundary submesh
   */
  void SetSubmeshField(const mfem::ParGridFunction& submesh_gridfn);

  /**
   * @brief Sets a new transfer object when the redecomp mesh has been updated
   * 
   * @param xfer Updated submesh to redecomp transfer object
   */
  void UpdateField(SubmeshRedecompTransfer& submesh_redecomp_xfer);

  /**
   * @brief Get the parent-linked boundary submesh grid function
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetSubmeshGridFn() const
  {
    return submesh_gridfn_;
  }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return mfem::GridFunction& 
   */
  mfem::GridFunction& GetRedecompGridFn() { return GetUpdateData().redecomp_gridfn_; }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompGridFn() const
  {
    return GetUpdateData().redecomp_gridfn_;
  }
  
  /**
   * @brief Get pointers to component arrays of the redecomp mesh grid function
   *
   * @return std::vector<const real*> of length 3
   *
   * @note Unused entries are nullptr.  Only the first entry is used with
   * frictionless contact.
   */
  std::vector<const real*> GetRedecompFieldPtrs() const;
  
  /**
   * @brief Get pointers to component arrays of a redecomp mesh grid function
   *
   * @param redecomp_gridfn Redecomp mesh grid function
   * @return std::vector<real*> of length 3
   *
   * @note Unused entries are nullptr.  Only the first entry is used with
   * frictionless contact.
   */
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_gridfn);

private:
  /**
   * @brief Creates and stores data that changes when the redecomp mesh is
   * updated
   */
  struct UpdateData
  {
    /**
     * @brief Construct a new UpdateData object
     *
     * @param submesh_redecomp_xfer Submesh to redecomp field transfer object
     * @param submesh_gridfn Grid function on the parent-linked boundary submesh
     */
    UpdateData(
      SubmeshRedecompTransfer& submesh_redecomp_xfer,
      const mfem::ParGridFunction& submesh_gridfn
    );
    
    /**
     * @brief Submesh to redecomp field transfer object
     */
    const SubmeshRedecompTransfer& submesh_redecomp_xfer_;

    /**
     * @brief Grid function values on the redecomp mesh
     */
    mfem::GridFunction redecomp_gridfn_;
  };
  
  /**
   * @brief Get the UpdateData object
   * 
   * @return UpdateData& 
   */
  UpdateData& GetUpdateData();

  /**
   * @brief Get the UpdateData object
   * 
   * @return const UpdateData& 
   */
  const UpdateData& GetUpdateData() const;

  /**
   * @brief Grid function on the parent-linked boundary submesh
   *
   * @note Stored as a reference wrapper so the reference can be updated
   */
  std::reference_wrapper<const mfem::ParGridFunction> submesh_gridfn_;

  /**
   * @brief UpdateData object created upon call to UpdateField()
   */
  std::unique_ptr<UpdateData> update_data_;
};

/**
 * @brief Stores MFEM and transfer data associated with parent vector fields
 * (displacement and velocity)
 */
class MfemMeshData
{
public:
  /**
   * @brief Construct a new MfemMeshData object
   *
   * @param mesh_id_1 Integer identifier for first Tribol registered mesh
   * @param mesh_id_2 Integer identifier for second Tribol registered mesh
   * @param parent_mesh Parent mesh, i.e. volume mesh of parent domain
   * @param current_coords Grid function on parent mesh holding current
   * coordinates
   * @param attributes_1 Mesh boundary attributes identifying surface elements
   * in the first Tribol registered mesh
   * @param attributes_2 Mesh boundary attributes identifying surface elements
   * in the second Tribol registered mesh
   */
  MfemMeshData(
    integer mesh_id_1,
    integer mesh_id_2,
    const mfem::ParMesh& parent_mesh,
    const mfem::ParGridFunction& current_coords,
    std::set<integer>&& attributes_1,
    std::set<integer>&& attributes_2
  );

  /**
   * @brief Get coordinate grid function on the parent mesh
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetParentCoords() const
  {
    return coords_.GetParentGridFn();
  }

  /**
   * @brief Sets a new coordinate grid function on the parent mesh
   * 
   * @param current_coords Coordinate grid function on the parent mesh
   */
  void SetParentCoords(const mfem::ParGridFunction& current_coords);

  /**
   * @brief Build a new redecomp mesh and update grid functions on the redecomp
   * mesh
   *
   * @note This method should be called after the coordinate grid function is
   * updated.
   */
  void UpdateMfemMeshData();

  /**
   * @brief Get the integer identifier for the first Tribol registered mesh
   *
   * @return integer 
   */
  integer GetMesh1ID() const { return mesh_id_1_; }

  /**
   * @brief Get the integer identifier for the second Tribol registered mesh
   *
   * @return integer 
   */
  integer GetMesh2ID() const { return mesh_id_2_; }

  /**
   * @brief Get the number of elements in the first Tribol registered mesh
   *
   * @return integer 
   */
  integer GetMesh1NE() const
  { 
    return GetUpdateData().conn_1_.size() / num_verts_per_elem_;
  }

  /**
   * @brief Get the number of elements in the second Tribol registered mesh
   *
   * @return integer 
   */
  integer GetMesh2NE() const
  {
    return GetUpdateData().conn_2_.size() / num_verts_per_elem_;
  }

  /**
   * @brief Get the total number of vertices in both Tribol registered meshes
   *
   * @return integer 
   */
  integer GetNV() const { return GetUpdateData().redecomp_mesh_.GetNV(); }

  /**
   * @brief Get the connectivity for the first Tribol registered mesh
   *
   * @return const IndexType* 
   */
  const IndexType* GetMesh1Conn() const
  {
    return GetUpdateData().conn_1_.data();
  }

  /**
   * @brief Get the connectivity for the second Tribol registered mesh
   *
   * @return const IndexType* 
   */
  const IndexType* GetMesh2Conn() const
  {
    return GetUpdateData().conn_2_.data();
  }

  /**
   * @brief Get the element type for both Tribol registered meshes
   *
   * @return integer 
   */
  integer GetElemType() const { return elem_type_; }

  /**
   * @brief Get pointers to component arrays of the coordinates on the
   * redecomp mesh
   *
   * @return std::vector<const real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  std::vector<const real*> GetRedecompCoordsPtrs() const 
  { 
    return coords_.GetRedecompFieldPtrs();
  }

  /**
   * @brief Get pointers to component arrays of the nodal response on the
   * redecomp mesh
   *
   * @return std::vector<real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  std::vector<real*> GetRedecompResponsePtrs()
  {
    return ParentField::GetRedecompFieldPtrs(redecomp_response_);
  }

  /**
   * @brief Get the nodal response grid function on the redecomp mesh
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompResponse() const
  {
    return redecomp_response_;
  }

  /**
   * @brief Get the nodal response vector on the parent mesh
   *
   * @param [out] r Pre-allocated, initialized mfem::Vector to which response
   * vector is added
   */
  void GetParentResponse(mfem::Vector& r) const;

  /**
   * @brief Get the parent to redecomp grid function transfer object
   * 
   * @return const ParentRedecompTransfer& 
   */
  const ParentRedecompTransfer& GetParentRedecompTransfer() const
  { 
    return GetUpdateData().vector_xfer_;
  }

  /**
   * @brief Add/replace the parent velocity grid function
   * 
   * @param velocity Velocity grid function on the parent mesh
   */
  void SetParentVelocity(const mfem::ParGridFunction& velocity);

  /**
   * @brief Determine if a velocity grid function has been set
   * 
   * @return true: Velocity grid function has been set
   * @return false: Velocity grid function has not been set
   */
  bool HasVelocity() const { return velocity_ != nullptr; }

  /**
   * @brief Get pointers to component arrays of the velocity on the RedecompMesh
   * 
   * @return std::vector<const real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  std::vector<const real*> GetRedecompVelocityPtrs() const
  {
    return velocity_->GetRedecompFieldPtrs();
  }

  /**
   * @brief Get the map from Tribol registered mesh 1 element indices to
   * redecomp mesh element indices
   *
   * @return const std::vector<integer>& 
   */
  const std::vector<integer>& GetElemMap1() const
  {
    return GetUpdateData().elem_map_1_;
  }

  /**
   * @brief Get the map from Tribol registered mesh 2 element indices to
   * redecomp mesh element indices
   *
   * @return const std::vector<integer>& 
   */
  const std::vector<integer>& GetElemMap2() const
  {
    return GetUpdateData().elem_map_2_;
  }

  /**
   * @brief Get the parent-linked boundary submesh containing both contact
   * surfaces
   *
   * @return mfem::ParSubMesh& 
   */
  mfem::ParSubMesh& GetSubmesh() { return submesh_; }

  /**
   * @brief Get the parent-linked boundary submesh containing both contact
   * surfaces
   *
   * @return const mfem::ParSubMesh& 
   */
  const mfem::ParSubMesh& GetSubmesh() const { return submesh_; }

  /**
   * @brief Get the LOR mesh containing both contact surfaces
   * 
   * @return mfem::ParMesh*
   *
   * @note nullptr if no refined mesh exists (polynomial order of parent is 1)
   */
  mfem::ParMesh* GetLORMesh() { return lor_mesh_.get(); }

  /**
   * @brief Get the LOR mesh containing both contact surfaces
   *
   * @return const mfem::ParMesh*
   *
   * @note nullptr if no refined mesh exists (polynomial order of parent is 1)
   */
  const mfem::ParMesh* GetLORMesh() const { return lor_mesh_.get(); }

  /**
   * @brief Get the redecomp mesh containing redecomposed contact surfaces
   * 
   * @return redecomp::RedecompMesh& 
   */
  redecomp::RedecompMesh& GetRedecompMesh()
  {
    return GetUpdateData().redecomp_mesh_;
  }

  /**
   * @brief Get the finite element space on the parent-linked boundary submesh
   *
   * @return const mfem::ParFiniteElementSpace& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_xfer_gridfn_.ParFESpace();
  }

  /**
   * @brief Get the finite element space on the LOR mesh
   *
   * @return const mfem::ParFiniteElementSpace* 
   *
   * @note nullptr if no LOR mesh exists (polynomial order of parent is 1)
   */
  const mfem::ParFiniteElementSpace* GetLORMeshFESpace() const
  {
    return submesh_lor_xfer_ ? submesh_lor_xfer_->GetLORGridFn().ParFESpace() : nullptr;
  }

  /**
   * @brief Set the number of element subdivisions per dimension on the LOR mesh
   *
   * @param lor_factor Number of element subdivisions per dimension
   */
  void SetLORFactor(integer lor_factor);
private:
  /**
   * @brief Creates and stores data that changes when the RedecompMesh is
   * updated
   */
  struct UpdateData
  {
    /**
     * @brief Construct a new UpdateData object
     *
     * @param submesh Parent-linked boundary submesh of contact elements
     * @param lor_mesh LOR mesh of contact elements (if using LOR; nullptr
     * otherwise)
     * @param parent_fes Vector finite element space on the original parent mesh
     * @param submesh_gridfn Grid function on the parent-linked boundary submesh
     * used to temporarily store variables being transferred
     * @param submesh_lor_xfer Submesh to LOR grid function transfer object (if
     * using LOR; nullptr otherwise)
     * @param attributes_1 Set of boundary attributes identifying elements in
     * the first Tribol registered mesh
     * @param attributes_2 Set of boundary attributes identifying elements in
     * the second Tribol registered mesh
     * @param num_verts_per_elem Number of vertices on each element
     */
    UpdateData(
      mfem::ParSubMesh& submesh,
      mfem::ParMesh* lor_mesh,
      const mfem::ParFiniteElementSpace& parent_fes,
      mfem::ParGridFunction& submesh_gridfn,
      SubmeshLORTransfer* submesh_lor_xfer,
      const std::set<integer>& attributes_1,
      const std::set<integer>& attributes_2,
      integer num_verts_per_elem
    );

    /**
     * @brief Redecomposed boundary element mesh
     */
    redecomp::RedecompMesh redecomp_mesh_;

    /**
     * @brief Parent mesh to redecomp mesh field transfer object
     */
    ParentRedecompTransfer vector_xfer_;

    /**
     * @brief Redecomp mesh element connectivity for the first Tribol registered
     * mesh
     */
    axom::Array<integer, 2> conn_1_;

    /**
     * @brief Redecomp mesh element connectivity for the second Tribol
     * registered mesh
     */
    axom::Array<integer, 2> conn_2_;

    /**
     * @brief Map from first Tribol registered mesh element indices to redecomp
     * mesh element indices
     */
    std::vector<integer> elem_map_1_;

    /**
     * @brief Map from second Tribol registered mesh element indices to redecomp
     * mesh element indices
     */
    std::vector<integer> elem_map_2_;

  private:
    /**
     * @brief Builds connectivity arrays and redecomp mesh to Tribol registered
     * mesh element maps
     *
     * @param attributes_1 Set of boundary attributes for the first Tribol
     * registered mesh
     * @param attributes_2 Set of boundary attributes for the second Tribol
     * registered mesh
     * @param num_verts_per_elem Number of vertices on each element
     */
    void UpdateConnectivity(
      const std::set<integer>& attributes_1,
      const std::set<integer>& attributes_2,
      integer num_verts_per_elem
    );
  };
  
  /**
   * @brief Get the UpdateData object
   * 
   * @return UpdateData& 
   */
  UpdateData& GetUpdateData();

  /**
   * @brief Get the UpdateData object
   * 
   * @return const UpdateData& 
   */
  const UpdateData& GetUpdateData() const;

  /**
   * @brief Create the parent-linked boundary submesh
   *
   * @param parent_mesh Parent mesh, i.e. volume mesh of parent domain
   * @param attributes_1 Mesh boundary attributes identifying surface elements
   * in the first Tribol registered mesh
   * @param attributes_2 Mesh boundary attributes identifying surface elements
   * in the second Tribol registered mesh
   * @return mfem::ParSubMesh 
   */
  static mfem::ParSubMesh CreateSubmesh(
    const mfem::ParMesh& parent_mesh,
    const std::set<integer>& attributes_1,
    const std::set<integer>& attributes_2
  );

  /**
   * @brief Create a grid function on the given parent-linked boundary submesh
   *
   * @param submesh Parent-linked boundary submesh
   * @param parent_fes Finite element space on the parent mesh
   * @return mfem::ParGridFunction 
   */
  static mfem::ParGridFunction CreateSubmeshGridFn(
    mfem::ParSubMesh& submesh,
    mfem::ParFiniteElementSpace& parent_fes
  );

  /**
   * @brief First mesh identifier
   */
  integer mesh_id_1_;

  /**
   * @brief Second mesh identifier
   */
  integer mesh_id_2_;

  /**
   * @brief Volume mesh of parent domain
   */
  const mfem::ParMesh& parent_mesh_;

  /**
   * @brief Mesh boundary attributes identifying first mesh
   */
  std::set<integer> attributes_1_;

  /**
   * @brief Mesh boundary attributes identifying second mesh
   */
  std::set<integer> attributes_2_;

  /**
   * @brief Submesh containing boundary elements of both contact surfaces
   */
  mfem::ParSubMesh submesh_;

  /**
   * @brief Coordinates grid function and transfer operators
   */
  ParentField coords_;

  /**
   * @brief Submesh grid function to temporarily hold values being transferred
   */
  mfem::ParGridFunction submesh_xfer_gridfn_;

  /**
   * @brief Refinement factor for refined mesh
   */
  integer lor_factor_;

  /**
   * @brief Contains LOR mesh if low-order refinement is being used; nullptr
   * otherwise
   */
  std::unique_ptr<mfem::ParMesh> lor_mesh_;

  /**
   * @brief Contains LOR mesh to submesh transfer operators if LOR is being
   * used; nullptr otherwise
   */
  std::unique_ptr<SubmeshLORTransfer> submesh_lor_xfer_;

  /**
   * @brief Type of elements on the contact meshes
   */
  InterfaceElementType elem_type_;

  /**
   * @brief Number of vertices on each element in the contact meshes
   */
  integer num_verts_per_elem_;

  /**
   * @brief Contains velocity grid function and transfer operators if set;
   * nullptr otherwise
   */
  std::unique_ptr<ParentField> velocity_;

  /**
   * @brief UpdateData object created upon call to UpdateMeshData()
   */
  std::unique_ptr<UpdateData> update_data_;

  /**
   * @brief Nodal response grid function on the redecomp mesh
   */
  mfem::GridFunction redecomp_response_;

  /**
   * @brief Merges two STL containers
   * 
   * @tparam T container type
   * @param container_1 First container
   * @param container_2 Second container
   * @return T merged container
   */
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

/**
 * @brief Stores MFEM and transfer data associated with parent-linked boundary
 * submesh pressure and gap fields 
 */
class MfemSubmeshData
{
public:
  /**
   * @brief Construct a new MfemSubmeshData object
   *
   * @param submesh Parent-linked boundary submesh
   * @param lor_mesh LOR mesh of contact surfaces (if using LOR; nullptr
   * otherwise)
   * @param pressure_fec Finite element collection of the pressure field
   * @param pressure_vdim Vector dimension of the pressure field
   */
  MfemSubmeshData(
    mfem::ParSubMesh& submesh,
    mfem::ParMesh* lor_mesh,
    std::unique_ptr<mfem::FiniteElementCollection> pressure_fec,
    integer pressure_vdim
  );

  /**
   * @brief Build a new transfer operator and update redecomp-level grid
   * functions
   *
   * @param redecomp_mesh Updated redecomp mesh
   */
  void UpdateMfemSubmeshData(redecomp::RedecompMesh& redecomp_mesh);

  /**
   * @brief Get pointers to component arrays of the pressure on the redecomp
   * mesh
   *
   * @return std::vector<const real*> of length 3
   *
   * @note Unused entries are nullptr.  Only the first entry is used with
   * frictionless contact.
   */
  std::vector<const real*> GetRedecompPressurePtrs() const
  {
    return pressure_.GetRedecompFieldPtrs();
  }

  /**
   * @brief Get the parent-linked boundary submesh pressure grid function
   * 
   * @return mfem::ParGridFunction& 
   */
  mfem::ParGridFunction& GetSubmeshPressure()
  {
    return submesh_pressure_;
  }

  /**
   * @brief Get the parent-linked boundary submesh pressure grid function
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetSubmeshPressure() const
  {
    return submesh_pressure_;
  }

  /**
   * @brief Get pointers to component arrays of the gap on the redecomp mesh
   * 
   * @return std::vector<real*> of length 3
   *
   * @note Unused entries are nullptr.  Only the first entry is used with
   * frictionless contact.
   */
  std::vector<real*> GetRedecompGapPtrs()
  {
    return PressureField::GetRedecompFieldPtrs(redecomp_gap_);
  }

  /**
   * @brief Get the gap grid function on the redecomp mesh
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompGap() const
  {
    return redecomp_gap_;
  }

  /**
   * @brief Get the gap vector on the parent-linked boundary submesh
   *
   * @param [out] g Un-initialized mfem::Vector holding the nodal gap values
   */
  void GetSubmeshGap(mfem::Vector& g) const;

  /**
   * @brief Get the parent-linked boundary submesh to redecomp mesh pressure
   * transfer object
   *
   * @return const SubmeshRedecompTransfer& 
   */
  const SubmeshRedecompTransfer& GetPressureTransfer() const
  {
    return GetUpdateData().pressure_xfer_;
  }

  /**
   * @brief Get the finite element space on the parent-linked boundary submesh
   *
   * @return const mfem::ParFiniteElementSpace& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_pressure_.ParFESpace();
  }

  /**
   * @brief Get the finite element space on the LOR mesh
   *
   * @return const mfem::ParFiniteElementSpace* or nullptr if no LOR mesh
   * (polynomial order of parent is 1)
   */
  const mfem::ParFiniteElementSpace* GetLORMeshFESpace() const
  {
    return submesh_lor_xfer_ ? submesh_lor_xfer_->GetLORGridFn().ParFESpace() : nullptr;
  }

private:
  /**
   * @brief Creates and stores data that changes when the redecomp mesh is
   * updated
   */
  struct UpdateData
  {
    /**
     * @brief Construct a new UpdateData object
     *
     * @param submesh_fes Pressure finite element space on the parent-linked
     * boundary submesh
     * @param submesh_lor_xfer Submesh to LOR grid function transfer object (if
     * using LOR; nullptr otherwise)
     * @param redecomp Redecomp mesh
     */
    UpdateData(
      mfem::ParFiniteElementSpace& submesh_fes,
      SubmeshLORTransfer* submesh_lor_xfer,
      redecomp::RedecompMesh& redecomp_mesh
    );

    /**
     * @brief Parent-linked boundary submesh to redecomp mesh field transfer
     * object
     */
    SubmeshRedecompTransfer pressure_xfer_;
  };

  /**
   * @brief Get the UpdateData object
   * 
   * @return UpdateData& 
   */
  UpdateData& GetUpdateData();

  /**
   * @brief Get the UpdateData object
   * 
   * @return const UpdateData& 
   */
  const UpdateData& GetUpdateData() const;

  /**
   * @brief Pressure grid function on the parent-linked boundary submesh
   */
  mfem::ParGridFunction submesh_pressure_;

  /**
   * @brief Pressure grid function and transfer operators
   */
  PressureField pressure_;

  /**
   * @brief Contains LOR mesh transfer operators if LOR is being used; nullptr
   * otherwise
   */
  std::unique_ptr<SubmeshLORTransfer> submesh_lor_xfer_;

  /**
   * @brief Gap grid function on the redecomp mesh
   */
  mfem::GridFunction redecomp_gap_;

  /**
   * @brief UpdateData object created upon call to UpdateSubmeshData()
   */
  std::unique_ptr<UpdateData> update_data_;
};

/**
 * @brief Simplifies transfer of Jacobian matrix data between MFEM and Tribol
 */
class MfemJacobianData
{
public:
  /**
   * @brief Construct a new MfemJacobianData object
   * 
   * @param parent_data MFEM data associated with displacement and velocity
   * @param submesh_data MFEM data associated with pressure and gap
   */
  MfemJacobianData(
    const MfemMeshData& parent_data,
    const MfemSubmeshData& submesh_data
  );

  /**
   * @brief Builds new transfer data after a new redecomp mesh has been built
   */
  void UpdateJacobianXfer();

  /**
   * @brief Returns Jacobian contributions as an mfem::BlockOperator
   * 
   * @param method_data Method data holding element Jacobians
   * @return std::unique_ptr<mfem::BlockOperator> 
   */
  std::unique_ptr<mfem::BlockOperator> GetMfemBlockJacobian(
    const MethodData& method_data
  ) const;

private:
  /**
   * @brief Creates and stores data that changes when the redecomp mesh is
   * updated
   */
  struct UpdateData
  {
    /**
     * @brief Construct a new UpdateData object
     * 
     * @param parent_data MFEM data associated with displacement and velocity
     * @param submesh_data MFEM data associated with pressure and gap
     */
    UpdateData(
      const MfemMeshData& parent_data,
      const MfemSubmeshData& submesh_data
    );

    /**
     * @brief Redecomp to parent-linked boundary submesh transfer operator
     */
    std::unique_ptr<redecomp::MatrixTransfer> submesh_redecomp_xfer_;
  };

  /**
   * @brief Get the UpdateData object
   * 
   * @return UpdateData& 
   */
  UpdateData& GetUpdateData();

  /**
   * @brief Get the UpdateData object
   * 
   * @return const UpdateData& 
   */
  const UpdateData& GetUpdateData() const;

  /**
   * @brief MFEM and transfer data associated with displacement and velocity
   */
  const MfemMeshData& parent_data_;

  /**
   * @brief MFEM and transfer data associated with pressure and gap
   */
  const MfemSubmeshData& submesh_data_;

  /**
   * @brief Array of offsets equal to number of displacement and pressure
   * degrees of freedom
   */
  mfem::Array<int> block_offsets_;

  /**
   * @brief List giving global parent vdof given the submesh vdof
   */
  mfem::Array<int> submesh2parent_vdof_list_;

  /**
   * @brief UpdateData object created upon calling UpdateMatrixXfer()
   */
  std::unique_ptr<UpdateData> update_data_;
};

} // end namespace tribol

#endif /* BUILD_REDECOMP */

#endif /* SRC_MESH_MFEMDATA_HPP_ */
