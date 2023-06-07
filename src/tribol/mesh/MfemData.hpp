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
#include "tribol/mesh/MethodCouplingData.hpp"

namespace tribol
{

/**
 * @brief Facilitates transfer to submesh low-order refined mesh
 *
 * This class simplifies transferring 1) a variable (e.g. displacement) from a
 * higher-order grid function on a coarse mesh to a low-order grid function on a
 * fine mesh and 2) a variable in the dual space (e.g. force) from a low-order
 * grid function on a fine mesh to a higher-order grid function on a coarse
 * mesh. The coarse mesh is intended to live on a mfem::ParSubMesh level, but
 * this is not required to use this class.
 */
class SubmeshLORTransfer
{
public:
  /**
   * @brief Construct a new SubmeshLORTransfer object
   * 
   * @param ho_submesh_fes Higher order finite element space on the coarse mesh
   * @param lor_submesh Refined mesh
   */
  SubmeshLORTransfer(
    mfem::ParFiniteElementSpace& ho_submesh_fes,
    mfem::ParMesh& lor_submesh
  );

  /**
   * @brief Transfers data from a higher-order grid function on a coarse mesh
   *
   * Data is transferred to the low-order grid function on a fine mesh that can
   * be accessed using GetLORGridFn().
   *
   * @param ho_src Source higher-order grid function on the coarse mesh
   */
  void TransferToLORGridFn(
    const mfem::ParGridFunction& ho_src
  );

  /**
   * @brief Transfers data to a higher-order grid function on a coarse mesh
   *
   * Data must be stored in the low-order grid function on a fine mesh accessed
   * using GetLORGridFn().
   *
   * @param ho_dst Destination higher-order grid function on the coarse mesh
   */
  void TransferFromLORGridFn(
    mfem::ParGridFunction& ho_dst
  ) const;

  /**
   * @brief Access the local low-order grid function on the refined mesh
   * 
   * @return mfem::ParGridFunction& 
   */
  mfem::ParGridFunction& GetLORGridFn() { return lor_gridfn_; }

  /**
   * @brief Access the local low-order grid function on the refined mesh
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetLORGridFn() const { return lor_gridfn_; }

private:
  /**
  * @brief Create low-order grid function on the refined mesh
  * 
  * @param lor_mesh Refined mesh
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
   * @brief Local low-order grid function on the refined mesh 
   */
  mfem::ParGridFunction lor_gridfn_;

  /**
   * @brief Low-order refined <-> higher-order coarse transfer object
   */
  mutable mfem::L2ProjectionGridTransfer lor_xfer_;
};

/**
 * @brief Facilitates transferring variables from the submesh to redecomp levels
 *
 * This class simplifies transferring field variables from/to a grid function on
 * an mfem::ParSubMesh to/from a grid function on a redecomp::RedecompMesh.  If
 * transferring to a low-order grid function on a refined mesh is also required,
 * this class will perform that transfer as well.
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
   * @param submesh_fes Finite element space on the (coarse) mfem::ParSubMesh
   * @param lor_xfer (Optional) low-order grid function, refined mesh transfer object
   * @param redecomp RCB redecomposed (refined) submesh
   */
  SubmeshRedecompTransfer(
    mfem::ParFiniteElementSpace& submesh_fes,
    SubmeshLORTransfer* lor_xfer,
    redecomp::RedecompMesh& redecomp
  );

  /**
   * @brief Transfer grid function on submesh to grid function on redecomp mesh
   * 
   * @param src Grid function on (coarse) submesh
   * @return mfem::GridFunction on (refined) redecomp mesh
   */
  mfem::GridFunction SubmeshToRedecomp(const mfem::ParGridFunction& src) const;

  /**
   * @brief Transfer grid function on redecomp mesh to grid function on submesh
   * 
   * @param src Grid function on (refined) redecomp mesh
   * @return mfem::ParGridFunction on (coarse) submesh
   */
  mfem::ParGridFunction RedecompToSubmesh(const mfem::GridFunction& src) const;

  /**
   * @brief Transfer grid function on redecomp mesh to grid function on submesh
   * 
   * @param src Grid function on (refined) redecomp mesh
   * @param dst Zero-valued grid function on (coarse) submesh
   */
  void RedecompToSubmesh(
    const mfem::GridFunction& src,
    mfem::ParGridFunction& dst
  ) const;

  /**
   * @brief Get the (coarse) Submesh associated with the SubmeshRedecompTransfer object
   * 
   * @return const mfem::ParSubMesh& 
   */
  const mfem::ParSubMesh& GetSubmesh() const 
  { 
    return static_cast<const mfem::ParSubMesh&>(*submesh_fes_.GetParMesh());
  }

  /**
   * @brief Set a new (refined) redecomp mesh associated with the existing (refined) submesh
   * 
   * @param redecomp RCB redecomposed (refined) submesh
   */
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);

private:
  /**
   * @brief Create a finite element space on the (refined) redecomp mesh
   * 
   * @param redecomp RCB redecomposed (refined) submesh
   * @param submesh_fes Original finite element space on the (refined) submesh
   * @return std::unique_ptr<mfem::FiniteElementSpace> 
   */
  static std::unique_ptr<mfem::FiniteElementSpace> CreateRedecompFESpace(
    redecomp::RedecompMesh& redecomp,
    mfem::ParFiniteElementSpace& submesh_fes
  );

  /**
   * @brief Finite element space on the (coarse) submesh 
   */
  mfem::ParFiniteElementSpace& submesh_fes_;

  /**
   * @brief Finite element space on the (refined) redecomp mesh
   */
  mutable std::unique_ptr<mfem::FiniteElementSpace> redecomp_fes_;

  /**
   * @brief Transfer object between low-order and higher-order grid functions
   */
  SubmeshLORTransfer* lor_xfer_;

  /**
   * @brief Transfer object between (refined) submesh and (refined) redecomp mesh
   */
  const redecomp::RedecompTransfer redecomp_xfer_;
};

/**
 * @brief Facilitates transferring variables from the original mesh to redecomp
 * levels
 *
 * This class simplifies transferring field variables from/to a grid function on
 * an mfem::ParMesh to/from a grid function on a redecomp::RedecompMesh.  If
 * transferring to a low-order grid function on a refined mesh is also required,
 * this class will perform that transfer as well.
 */
class ParentRedecompTransfer
{
public:
  /**
   * @brief Construct a new ParentRedecompTransfer object
   *
   * @param parent_fes Finite element space on the original (coarse) mesh
   * @param submesh_gridfn Grid function on the (coarse) submesh used to
   * temporarily store variables being transferred
   * @param lor_xfer (Optional) Low-order grid function, refined mesh transfer
   * object
   * @param redecomp RCB redecomposed (refined) submesh
   */
  ParentRedecompTransfer(
    const mfem::ParFiniteElementSpace& parent_fes,
    mfem::ParGridFunction& submesh_gridfn,
    SubmeshLORTransfer* lor_xfer,
    redecomp::RedecompMesh& redecomp
  );

  /**
   * @brief Transfer grid function on original mesh to grid function on redecomp
   * mesh
   *
   * @param src Grid function on (coarse) original mesh
   * @return mfem::GridFunction on (refined) redecomp mesh
   */
  mfem::GridFunction ParentToRedecomp(const mfem::ParGridFunction& src) const;
  
  /**
   * @brief Transfer grid function on redecomp mesh to grid function on original mesh
   * 
   * @param src Grid function on (refined) redecomp mesh
   * @return mfem::ParGridFunction on (coarse) original mesh
   */
  mfem::ParGridFunction RedecompToParent(const mfem::GridFunction& src) const;

  /**
   * @brief Set a new (refined) redecomp mesh associated with the existing (refined) submesh
   * 
   * @param redecomp RCB redecomposed (refined) submesh
   */
  void UpdateRedecomp(redecomp::RedecompMesh& redecomp);

  /**
   * @brief Get the (coarse) Submesh finite element space associated with this
   * transfer object
   *
   * @return const mfem::ParSubMesh& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_gridfn_.ParFESpace();
  }

private:
  /**
   * @brief Finite element space on the (coarse) original mesh
   */
  mutable mfem::ParFiniteElementSpace parent_fes_;

  /**
   * @brief Grid function on the (coarse) submesh
   */
  mfem::ParGridFunction& submesh_gridfn_;

  /**
   * @brief Object to transfer variables to/from the (coarse) submesh level
   * from/to the (refined) redecomp level
   */
  SubmeshRedecompTransfer redecomp_xfer_;
};

/**
 * @brief Vector field variable that lives on the original, parent mesh
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
   * @param parent Grid function on the original (parent) mesh
   */
  ParentField(const mfem::ParGridFunction& parent);
  
  /**
   * @brief Set a new grid function on the original (parent) mesh
   * 
   * @param parent Grid function on the original (parent) mesh
   */
  void SetParentField(const mfem::ParGridFunction& parent);

  /**
   * @brief Set a new transfer object when the redecomp mesh has been updated
   *
   * @param xfer Updated parent to redecomp transfer object
   */
  void UpdateField(const ParentRedecompTransfer& xfer);

  /**
   * @brief Get the parent grid function
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetParentField() const
  {
    return parent_;
  }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return mfem::GridFunction& 
   */
  mfem::GridFunction& GetRedecompField() { return GetUpdateData().redecomp_; }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompField() const
  {
    return GetUpdateData().redecomp_;
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
   * @param redecomp_fn Redecomp mesh grid function
   * @return std::vector<real*> of length 3
   *
   * @note The third entry is nullptr in two dimensions
   */
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn);

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
     * @param xfer Parent to redecomp field transfer object
     * @param parent Grid function on the original, parent mesh
     */
    UpdateData(
      const ParentRedecompTransfer& xfer,
      const mfem::ParGridFunction& parent
    );

    /**
     * @brief Parent to redecomp field transfer object
     */
    const ParentRedecompTransfer& xfer_;

    /**
     * @brief Grid function values on the redecomp mesh
     */
    mfem::GridFunction redecomp_;
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
   * @brief Grid function on the parent, original mesh
   *
   * @note Stored as a reference wrapper so the reference can be updated
   */
  std::reference_wrapper<const mfem::ParGridFunction> parent_;

  /**
   * @brief UpdateData object created upon call to UpdateField()
   */
  std::unique_ptr<UpdateData> update_data_;
};

/**
 * @brief Stores a pressure field variable that lives on the (coarse) submesh
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
   * @param submesh Grid function on the (coarse) submesh
   */
  PressureField(const mfem::ParGridFunction& submesh);

  /**
   * @brief Sets a new grid function on the (coarse) submesh
   * 
   * @param submesh Grid function on the (coarse) submesh
   */
  void SetSubmeshField(const mfem::ParGridFunction& submesh);

  /**
   * @brief Sets a new transfer object when the redecomp mesh has been updated
   * 
   * @param xfer Updated submesh to redecomp transfer object
   */
  void UpdateField(const SubmeshRedecompTransfer& xfer);

  /**
   * @brief Get the submesh grid function
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetSubmeshField() const
  {
    return submesh_;
  }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return mfem::GridFunction& 
   */
  mfem::GridFunction& GetRedecompField() { return GetUpdateData().redecomp_; }

  /**
   * @brief Get the redecomp mesh grid function
   * 
   * @return const mfem::GridFunction& 
   */
  const mfem::GridFunction& GetRedecompField() const
  {
    return GetUpdateData().redecomp_;
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
   * @brief Get pointers to component arrays of the redecomp mesh grid function
   *
   * @param redecomp_fn Redecomp mesh grid function
   * @return std::vector<real*> of length 3
   *
   * @note Unused entries are nullptr.  Only the first entry is used with
   * frictionless contact.
   */
  static std::vector<real*> GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn);

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
     * @param xfer Submesh to redecomp field transfer object
     * @param submesh Grid function on the (coarse) submesh
     */
    UpdateData(
      const SubmeshRedecompTransfer& xfer,
      const mfem::ParGridFunction& submesh
    );
    
    /**
     * @brief Submesh to redecomp field transfer object
     */
    const SubmeshRedecompTransfer& xfer_;

    /**
     * @brief Grid function values on the redecomp mesh
     */
    mfem::GridFunction redecomp_;
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
   * @brief Grid function on the submesh
   *
   * @note Stored as a reference wrapper so the reference can be updated
   */
  std::reference_wrapper<const mfem::ParGridFunction> submesh_;

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
   * @param mesh Volume mesh of parent domain
   * @param current_coords Grid function (on mesh) holding current coordinates
   * @param attributes_1 Mesh boundary attributes identifying first mesh
   * @param attributes_2 Mesh boundary attributes identifying second mesh
   */
  MfemMeshData(
    integer mesh_id_1,
    integer mesh_id_2,
    mfem::ParMesh& mesh,
    const mfem::ParGridFunction& current_coords,
    const std::set<integer>& attributes_1,
    const std::set<integer>& attributes_2
  );

  /**
   * @brief Get coordinate grid function on the parent mesh
   * 
   * @return const mfem::ParGridFunction& 
   */
  const mfem::ParGridFunction& GetParentCoords() const
  {
    return coords_.GetParentField();
  }

  /**
   * @brief Sets a new coordinate grid function on the parent mesh
   * 
   * @param current_coords Coordinate grid function on the parent mesh
   */
  void SetParentCoords(const mfem::ParGridFunction& current_coords);

  /**
   * @brief Build a new redecomp mesh and update grid functions
   *
   * @note This method should be called after the coordinate grid function has
   * changed.
   */
  void UpdateMeshData();

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
   * @brief Get the number of elements on the first mesh
   * 
   * @return integer 
   */
  integer GetMesh1NE() const
  { 
    return GetUpdateData().conn_1_.size() / num_verts_per_elem_;
  }

  /**
   * @brief Get the number of elements on the second mesh
   * 
   * @return integer 
   */
  integer GetMesh2NE() const
  {
    return GetUpdateData().conn_2_.size() / num_verts_per_elem_;
  }

  /**
   * @brief Get the total number of vertices in both meshes
   * 
   * @return integer 
   */
  integer GetNV() const { return GetUpdateData().redecomp_.GetNV(); }

  /**
   * @brief Get the connectivity for the first mesh
   * 
   * @return const IndexType* 
   */
  const IndexType* GetMesh1Conn() const
  {
    return GetUpdateData().conn_1_.data();
  }

  /**
   * @brief Get the connectivity for the second mesh
   * 
   * @return const IndexType* 
   */
  const IndexType* GetMesh2Conn() const
  {
    return GetUpdateData().conn_2_.data();
  }

  /**
   * @brief Get the element type for the meshes
   * 
   * @return integer 
   */
  integer GetElemType() const { return elem_type_; }

  /**
   * @brief Get pointers to component arrays of the coordinates on the redecomp mesh
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
   * @brief Get pointers to component arrays of the nodal response on the redecomp mesh
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
   * @brief Get the nodal response grid function on the parent mesh
   * 
   * @return mfem::ParGridFunction 
   */
  mfem::ParGridFunction GetParentResponse() const
  {
    return GetParentTransfer().RedecompToParent(redecomp_response_);
  }

  /**
   * @brief Get the parent to redecomp transfer object
   * 
   * @return const ParentRedecompTransfer& 
   */
  const ParentRedecompTransfer& GetParentTransfer() const
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
   * @brief Get pointers to component arrays of the velocity on the redecomp mesh
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
   * @brief Get the map from mesh 1 element indices to redecomp mesh element
   * indices
   *
   * @return const std::vector<integer>& 
   */
  const std::vector<integer>& GetElemMap1() const
  {
    return GetUpdateData().elem_map_1_;
  }

  /**
   * @brief Get the map from mesh 2 element indices to redecomp mesh element
   * indices
   *
   * @return const std::vector<integer>& 
   */
  const std::vector<integer>& GetElemMap2() const
  {
    return GetUpdateData().elem_map_2_;
  }

  /**
   * @brief Get the (coarse) boundary submesh containing contact surfaces
   * 
   * @return mfem::ParSubMesh& 
   */
  mfem::ParSubMesh& GetSubmesh() { return submesh_; }

  /**
   * @brief Get the (coarse) boundary submesh containing contact surfaces
   * 
   * @return const mfem::ParSubMesh& 
   */
  const mfem::ParSubMesh& GetSubmesh() const { return submesh_; }

  /**
   * @brief Get the (refined) boundary submesh containing contact surfaces
   * 
   * @return mfem::ParMesh*
   */
  mfem::ParMesh* GetLORSubmesh() { return lor_submesh_.get(); }

  /**
   * @brief Get the refined boundary submesh containing contact surfaces
   * 
   * @return const mfem::ParMesh*
   *
   * @note nullptr if no refined mesh exists
   */
  const mfem::ParMesh* GetLORSubmesh() const { return lor_submesh_.get(); }

  /**
   * @brief Get the (refined) redecomposed boundary submesh
   * 
   * @return redecomp::RedecompMesh& 
   */
  redecomp::RedecompMesh& GetRedecompMesh()
  {
    return GetUpdateData().redecomp_;
  }

  /**
   * @brief Get the finite element space on the (coarse) boundary submesh
   * 
   * @return const mfem::ParFiniteElementSpace& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_xfer_gridfn_.ParFESpace();
  }

  /**
   * @brief Get the finite element space on the refined boundary submesh
   * 
   * @return const mfem::ParFiniteElementSpace* 
   *
   * @note nullptr if no refined mesh exists
   */
  const mfem::ParFiniteElementSpace* GetLORSubmeshFESpace() const
  {
    return lor_xfer_ ? lor_xfer_->GetLORGridFn().ParFESpace() : nullptr;
  }

  /**
   * @brief Set the number of element subdivisions per dimension on the refined
   * mesh
   *
   * @param lor_factor Number of element subdivisions per dimension
   */
  void SetLowOrderRefinedFactor(integer lor_factor);
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
     * @param submesh Coarse boundary submesh of contact elements
     * @param lor_submesh (Optional) Refined boundary submesh of contact elements
     * @param parent_fes Vector finite element space on the original parent mesh
     * @param submesh_gridfn Grid function on the (coarse) submesh used to
     * temporarily store variables being transferred
     * @param lor_xfer (Optional) Low-order grid function, refined mesh transfer
     * object
     * @param attributes_1 Set of boundary attributes for the first mesh
     * @param attributes_2 Set of boundary attributes for the second mesh
     * @param num_verts_per_elem Number of vertices on each element
     */
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

    /**
     * @brief Redecomposed (refined) boundary element mesh
     */
    redecomp::RedecompMesh redecomp_;

    /**
     * @brief Parent to redecomp field transfer object
     */
    ParentRedecompTransfer vector_xfer_;

    /**
     * @brief Redecomp mesh element connectivity for the first mesh
     */
    axom::Array<integer, 2> conn_1_;

    /**
     * @brief Redecomp mesh element connectivity for the second mesh
     */
    axom::Array<integer, 2> conn_2_;

    /**
     * @brief Map from mesh 1 element indicies to redecomp mesh element indices
     */
    std::vector<integer> elem_map_1_;

    /**
     * @brief Map from mesh 2 element indicies to redecomp mesh element indices
     */
    std::vector<integer> elem_map_2_;

  private:
    /**
     * @brief Builds connectivity arrays and redecomp to tribol element maps
     * 
     * @param attributes_1 Set of boundary attributes for the first mesh
     * @param attributes_2 Set of boundary attributes for the second mesh
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
   * @brief Create the boundary element submesh
   * 
   * @param parent Volume mesh of parent domain
   * @param attributes_1 Mesh boundary attributes identifying first mesh
   * @param attributes_2 Mesh boundary attributes identifying second mesh
   * @return mfem::ParSubMesh 
   */
  static mfem::ParSubMesh CreateSubmesh(
    const mfem::ParMesh& parent,
    const std::set<integer>& attributes_1,
    const std::set<integer>& attributes_2
  );

  /**
   * @brief Create a grid function on the given submesh
   * 
   * @param submesh (Coarse) submesh of boundary elements on the contact surface
   * @param parent_fes Finite element space on the parent (volume) mesh
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
  mfem::ParMesh& mesh_;

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
   * @brief (Coarse) submesh grid function to temporarily hold values being
   * transferred
   */
  mfem::ParGridFunction submesh_xfer_gridfn_;

  /**
   * @brief Refinement factor for refined mesh
   */
  integer lor_factor_;

  /**
   * @brief Contains refined submesh if low-order refinement is being used;
   * nullptr otherwise
   */
  std::unique_ptr<mfem::ParMesh> lor_submesh_;

  /**
   * @brief Contains low-order refined submesh transfer operators if low-order
   * refinement is being used; nullptr otherwise
   */
  std::unique_ptr<SubmeshLORTransfer> lor_xfer_;

  /**
   * @brief Type of elements on the contact mesh
   */
  InterfaceElementType elem_type_;

  /**
   * @brief Number of vertices on each element in the contact mesh
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
 * @brief Stores MFEM and transfer data associated with submesh pressure and gap
 * fields 
 */
class MfemSubmeshData
{
public:
  /**
   * @brief Construct a new MfemSubmeshData object
   *
   * @param submesh (Coarse) boundary element mesh containing contact surfaces
   * @param lor_submesh (Optional) refined boundary element mesh containing
   * contact surfaces
   * @param pressure_fec Finite element collection of the pressure field
   * @param pressure_vdim Vector dimension of the pressure field
   */
  MfemSubmeshData(
    mfem::ParSubMesh& submesh,
    mfem::ParMesh* lor_submesh,
    std::unique_ptr<mfem::FiniteElementCollection> pressure_fec,
    integer pressure_vdim
  );

  /**
   * @brief Build a new transfer operator and update grid functions
   * 
   * @param redecomp Updated (refined) redecomp mesh
   */
  void UpdateSubmeshData(redecomp::RedecompMesh& redecomp);

  /**
   * @brief Get pointers to component arrays of the pressure on the redecomp mesh
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
   * @brief Get the (coarse) submesh pressure grid function
   * 
   * @return mfem::ParGridFunction& 
   */
  mfem::ParGridFunction& GetSubmeshPressure()
  {
    return submesh_pressure_;
  }

  /**
   * @brief Get the (coarse) submesh pressure grid function
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
   * @brief Get the gap grid function on the (coarse) submesh
   * 
   * @return mfem::ParGridFunction 
   */
  mfem::ParGridFunction GetSubmeshGap() const
  {
    return GetPressureTransfer().RedecompToSubmesh(redecomp_gap_);
  }

  /**
   * @brief Get the submesh to redecomp pressure transfer object
   * 
   * @return const SubmeshRedecompTransfer& 
   */
  const SubmeshRedecompTransfer& GetPressureTransfer() const
  {
    return GetUpdateData().pressure_xfer_;
  }

  /**
   * @brief Get the finite element space on the (coarse) submesh
   * 
   * @return const mfem::ParFiniteElementSpace& 
   */
  const mfem::ParFiniteElementSpace& GetSubmeshFESpace() const
  {
    return *submesh_pressure_.ParFESpace();
  }

  /**
   * @brief Get the finite element space on the refined submesh
   *
   * @return const mfem::ParFiniteElementSpace* or nullptr if no refined submesh
   */
  const mfem::ParFiniteElementSpace* GetLORSubmeshFESpace() const
  {
    return lor_xfer_ ? lor_xfer_->GetLORGridFn().ParFESpace() : nullptr;
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
     * @param submesh_fes Pressure finite element space on the (coarse) submesh
     * @param lor_xfer (Optional) Low-order grid function, refined mesh transfer
     * object
     * @param redecomp Updated (refined) redecomp mesh
     */
    UpdateData(
      mfem::ParFiniteElementSpace& submesh_fes,
      SubmeshLORTransfer* lor_xfer,
      redecomp::RedecompMesh& redecomp
    );

    /**
     * @brief Submesh to redecomp field transfer object
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
   * @brief Pressure grid function on the (coarse) submesh
   */
  mfem::ParGridFunction submesh_pressure_;

  /**
   * @brief Pressure grid function and transfer operators
   */
  PressureField pressure_;

  /**
   * @brief Contains low-order refined submesh transfer operators if low-order
   * refinement is being used; nullptr otherwise
   */
  std::unique_ptr<SubmeshLORTransfer> lor_xfer_;

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
 * @brief Simplifies transfer of matrix data between MFEM and Tribol
 */
class MfemMatrixData
{
public:
  /**
   * @brief Construct a new MfemMatrixData object
   * 
   * @param parent_data MFEM data associated with displacement and velocity
   * @param submesh_data MFEM data associated with pressure and gap
   */
  MfemMatrixData(
    const MfemMeshData& parent_data,
    const MfemSubmeshData& submesh_data
  );

  /**
   * @brief Builds new transfer data after a new redecomp mesh has been built
   */
  void UpdateMatrixXfer();

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
     * @brief Redecomp to (refined) submesh transfer operator
     */
    std::unique_ptr<redecomp::MatrixTransfer> matrix_xfer_;
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

#endif /* SRC_MESH_MFEMDATA_HPP_ */
