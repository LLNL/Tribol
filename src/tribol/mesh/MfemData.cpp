// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)


#include "tribol/mesh/MfemData.hpp"

#ifdef BUILD_REDECOMP

#include "axom/slic.hpp"

namespace tribol
{

SubmeshLORTransfer::SubmeshLORTransfer(
  mfem::ParFiniteElementSpace& submesh_fes,
  mfem::ParMesh& lor_mesh
)
: lor_gridfn_ { CreateLORGridFunction(
    lor_mesh,
    std::make_unique<mfem::H1_FECollection>(1, lor_mesh.SpaceDimension()),
    submesh_fes.GetVDim()
  ) },
  lor_xfer_ { submesh_fes, *lor_gridfn_.ParFESpace() }
{
  SLIC_WARNING_ROOT("LOR support is experimental at this time.");
}

void SubmeshLORTransfer::TransferToLORGridFn(
  const mfem::ParGridFunction& submesh_src
)
{
  lor_xfer_.ForwardOperator().Mult(submesh_src, lor_gridfn_);
}

void SubmeshLORTransfer::TransferFromLORGridFn(
  mfem::ParGridFunction& submesh_dst
) const
{
  lor_xfer_.ForwardOperator().MultTranspose(lor_gridfn_, submesh_dst);
}

mfem::ParGridFunction SubmeshLORTransfer::CreateLORGridFunction(
  mfem::ParMesh& lor_mesh,
  std::unique_ptr<mfem::FiniteElementCollection> lor_fec,
  integer vdim
)
{
  mfem::ParGridFunction lor_gridfn {
    new mfem::ParFiniteElementSpace(
      &lor_mesh,
      lor_fec.get(),
      vdim,
      mfem::Ordering::byNODES
    )
  };
  lor_gridfn.MakeOwner(lor_fec.release());
  return lor_gridfn;
}

SubmeshRedecompTransfer::SubmeshRedecompTransfer(
  mfem::ParFiniteElementSpace& submesh_fes,
  SubmeshLORTransfer* submesh_lor_xfer,
  redecomp::RedecompMesh& redecomp_mesh
)
: submesh_fes_ { submesh_fes },
  redecomp_fes_ { submesh_lor_xfer ?
    CreateRedecompFESpace(redecomp_mesh, *submesh_lor_xfer->GetLORGridFn().ParFESpace()) :
    CreateRedecompFESpace(redecomp_mesh, submesh_fes_) },
  submesh_lor_xfer_ { submesh_lor_xfer },
  redecomp_xfer_ { } // default (element transfer) constructor
{
  // make sure submesh_fes is a submesh and redecomp's parent is submesh_fes's
  // submesh
  SLIC_ERROR_ROOT_IF(
    !mfem::ParSubMesh::IsParSubMesh(submesh_fes_.GetParMesh()),
    "submesh_fes must be on a ParSubMesh."
  );
  SLIC_ERROR_ROOT_IF(
    !submesh_lor_xfer &&
      &redecomp_mesh.getParent() != submesh_fes_.GetParMesh(),
    "redecomp's parent must match the submesh_fes ParMesh."
  );
  SLIC_ERROR_ROOT_IF(
    submesh_lor_xfer &&
      &redecomp_mesh.getParent() != submesh_lor_xfer->GetLORGridFn().ParFESpace()->GetParMesh(),
    "redecomp's parent must match the submesh_fes ParMesh."
  );
}

void SubmeshRedecompTransfer::SubmeshToRedecomp(
  const mfem::ParGridFunction& submesh_src,
  mfem::GridFunction& redecomp_dst
) const
{
  auto src_ptr = &submesh_src;
  if (submesh_lor_xfer_)
  {
    submesh_lor_xfer_->GetLORGridFn() = 0.0;
    submesh_lor_xfer_->TransferToLORGridFn(submesh_src);
    src_ptr = &submesh_lor_xfer_->GetLORGridFn();
  }
  redecomp_xfer_.TransferToSerial(*src_ptr, redecomp_dst);
}

void SubmeshRedecompTransfer::RedecompToSubmesh(
  const mfem::GridFunction& redecomp_src,
  mfem::ParGridFunction& submesh_dst
) const
{
  auto dst_ptr = &submesh_dst;
  if (submesh_lor_xfer_)
  {
    submesh_lor_xfer_->GetLORGridFn() = 0.0;
    dst_ptr = &submesh_lor_xfer_->GetLORGridFn();
  }
  redecomp_xfer_.TransferToParallel(redecomp_src, *dst_ptr);
  if (submesh_lor_xfer_)
  {
    submesh_lor_xfer_->TransferFromLORGridFn(submesh_dst);
  }
}

std::unique_ptr<mfem::FiniteElementSpace> SubmeshRedecompTransfer::CreateRedecompFESpace(
  redecomp::RedecompMesh& redecomp_mesh,
  mfem::ParFiniteElementSpace& submesh_fes
)
{
  return std::make_unique<mfem::FiniteElementSpace>(
    &redecomp_mesh,
    submesh_fes.FEColl(),
    submesh_fes.GetVDim(),
    mfem::Ordering::byNODES
  );
}

ParentRedecompTransfer::ParentRedecompTransfer(
  const mfem::ParFiniteElementSpace& parent_fes,
  mfem::ParGridFunction& submesh_gridfn,
  SubmeshLORTransfer* submesh_lor_xfer,
  redecomp::RedecompMesh& redecomp_mesh
)
: parent_fes_ { parent_fes },
  submesh_gridfn_ { submesh_gridfn },
  submesh_redecomp_xfer_ { 
    *submesh_gridfn_.ParFESpace(),
    submesh_lor_xfer, 
    redecomp_mesh
  }
{
  // Note: this is checked in the SubmeshRedecompTransfer constructor
  // SLIC_ERROR_ROOT_IF(
  //   !mfem::ParSubMesh::IsParSubMesh(submesh_gridfn_.ParFESpace()->GetParMesh()),
  //   "submesh_gridfn_ must be associated with an mfem::ParSubMesh."
  // );
  SLIC_ERROR_ROOT_IF(
    submesh_redecomp_xfer_.GetSubmesh().GetParent() != parent_fes_.GetParMesh(),
    "submesh_gridfn's parent mesh must match the parent_fes ParMesh."
  );
}

void ParentRedecompTransfer::ParentToRedecomp(
  const mfem::ParGridFunction& parent_src,
  mfem::GridFunction& redecomp_dst
) const
{
  submesh_gridfn_ = 0.0;
  submesh_redecomp_xfer_.GetSubmesh().Transfer(parent_src, submesh_gridfn_);
  submesh_redecomp_xfer_.SubmeshToRedecomp(submesh_gridfn_, redecomp_dst);
}

void ParentRedecompTransfer::RedecompToParent(
  const mfem::GridFunction& redecomp_src,
  mfem::ParGridFunction& parent_dst
) const
{
  submesh_gridfn_ = 0.0;
  submesh_redecomp_xfer_.RedecompToSubmesh(redecomp_src, submesh_gridfn_);
  submesh_redecomp_xfer_.GetSubmesh().Transfer(submesh_gridfn_, parent_dst);
}

ParentField::ParentField(
  const mfem::ParGridFunction& parent_gridfn
)
: parent_gridfn_ { parent_gridfn }
{}

void ParentField::SetParentGridFn(const mfem::ParGridFunction& parent_gridfn)
{
  parent_gridfn_ = parent_gridfn;
  update_data_.reset(nullptr);
}

void ParentField::UpdateField(ParentRedecompTransfer& parent_redecomp_xfer)
{
  update_data_ = std::make_unique<UpdateData>(parent_redecomp_xfer, parent_gridfn_);
}

std::vector<const real*> ParentField::GetRedecompFieldPtrs() const
{
  auto data_ptrs = std::vector<const real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(GetRedecompGridFn().FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &GetRedecompGridFn()(GetRedecompGridFn().FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

std::vector<real*> ParentField::GetRedecompFieldPtrs(mfem::GridFunction& redecomp_gridfn)
{
  auto data_ptrs = std::vector<real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(redecomp_gridfn.FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &redecomp_gridfn(redecomp_gridfn.FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

ParentField::UpdateData& ParentField::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const ParentField::UpdateData& ParentField::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

ParentField::UpdateData::UpdateData(
  ParentRedecompTransfer& parent_redecomp_xfer,
  const mfem::ParGridFunction& parent_gridfn
)
: parent_redecomp_xfer_ { parent_redecomp_xfer },
  redecomp_gridfn_ { &parent_redecomp_xfer.GetRedecompFESpace() }
{
  redecomp_gridfn_ = 0.0;
  parent_redecomp_xfer_.ParentToRedecomp(parent_gridfn, redecomp_gridfn_);
}

PressureField::PressureField(
  const mfem::ParGridFunction& submesh_gridfn
)
: submesh_gridfn_ { submesh_gridfn }
{}

void PressureField::SetSubmeshField(const mfem::ParGridFunction& submesh_gridfn)
{
  submesh_gridfn_ = submesh_gridfn;
  update_data_.reset(nullptr);
}

void PressureField::UpdateField(SubmeshRedecompTransfer& submesh_redecomp_xfer)
{
  update_data_ = std::make_unique<UpdateData>(submesh_redecomp_xfer, submesh_gridfn_);
}

std::vector<const real*> PressureField::GetRedecompFieldPtrs() const
{
  auto data_ptrs = std::vector<const real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(GetRedecompGridFn().FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &GetRedecompGridFn()(GetRedecompGridFn().FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

std::vector<real*> PressureField::GetRedecompFieldPtrs(mfem::GridFunction& redecomp_gridfn)
{
  auto data_ptrs = std::vector<real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(redecomp_gridfn.FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &redecomp_gridfn(redecomp_gridfn.FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

PressureField::UpdateData& PressureField::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const PressureField::UpdateData& PressureField::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

PressureField::UpdateData::UpdateData(
  SubmeshRedecompTransfer& submesh_redecomp_xfer,
  const mfem::ParGridFunction& submesh_gridfn
)
: submesh_redecomp_xfer_ { submesh_redecomp_xfer },
  redecomp_gridfn_ { &submesh_redecomp_xfer.GetRedecompFESpace() }
{
  redecomp_gridfn_ = 0.0;
  submesh_redecomp_xfer_.SubmeshToRedecomp(submesh_gridfn, redecomp_gridfn_);
}

MfemMeshData::MfemMeshData(
  integer mesh_id_1,
  integer mesh_id_2,
  const mfem::ParMesh& parent_mesh,
  const mfem::ParGridFunction& current_coords,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2
)
: mesh_id_1_ { mesh_id_1 },
  mesh_id_2_ { mesh_id_2 },
  parent_mesh_ { parent_mesh },
  attributes_1_ { attributes_1 },
  attributes_2_ { attributes_2 },
  submesh_ { CreateSubmesh(parent_mesh_, attributes_1_, attributes_2_) },
  coords_ { current_coords },
  lor_factor_ { 0 }
{
  // create submesh grid function
  std::unique_ptr<mfem::FiniteElementCollection> submesh_fec {
    current_coords.ParFESpace()->FEColl()->Clone(
      current_coords.ParFESpace()->FEColl()->GetOrder()
    )
  };
  submesh_xfer_gridfn_.SetSpace(
    new mfem::ParFiniteElementSpace(
      &submesh_,
      submesh_fec.get(),
      current_coords.ParFESpace()->GetVDim(),
      mfem::Ordering::byNODES
    )
  );
  submesh_xfer_gridfn_.MakeOwner(submesh_fec.release());

  // build LOR submesh
  if (current_coords.FESpace()->FEColl()->GetOrder() > 1)
  {
    SetLORFactor(current_coords.FESpace()->FEColl()->GetOrder());
  }

  // set the element type
  mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  if (submesh_.GetNE() > 0)
  {
    element_type = submesh_.GetElementType(0);
  }

  switch (element_type) 
  {
    case mfem::Element::SEGMENT:
      elem_type_ = LINEAR_EDGE;
      break;
    case mfem::Element::TRIANGLE:
      elem_type_ = LINEAR_TRIANGLE;
      break;
    case mfem::Element::QUADRILATERAL:
      elem_type_ = LINEAR_QUAD;
      break;
    case mfem::Element::TETRAHEDRON:
      elem_type_ = LINEAR_TET;
      break;
    case mfem::Element::HEXAHEDRON:
      elem_type_ = LINEAR_HEX;
      break;

    case mfem::Element::POINT:
      SLIC_ERROR_ROOT("Unsupported element type!");
      break;

    default:
      SLIC_ERROR_ROOT("Unknown element type!");
      break;
  }

  num_verts_per_elem_ = mfem::Geometry::NumVerts[element_type];
}

void MfemMeshData::SetParentCoords(const mfem::ParGridFunction& current_coords)
{
  coords_.SetParentGridFn(current_coords);
}

void MfemMeshData::UpdateMfemMeshData()
{
  update_data_ = std::make_unique<UpdateData>(
    submesh_,
    lor_mesh_.get(),
    *coords_.GetParentGridFn().ParFESpace(),
    submesh_xfer_gridfn_,
    submesh_lor_xfer_.get(),
    attributes_1_, 
    attributes_2_, 
    num_verts_per_elem_
  );
  coords_.UpdateField(update_data_->vector_xfer_);
  redecomp_response_.SetSpace(coords_.GetRedecompGridFn().FESpace());
  redecomp_response_ = 0.0;
  if (velocity_)
  {
    velocity_->UpdateField(update_data_->vector_xfer_);
  }
}

void MfemMeshData::GetParentResponse(mfem::Vector& r) const
{
    mfem::ParGridFunction r_gridfn(coords_.GetParentGridFn().ParFESpace(), r);
    GetParentRedecompTransfer().RedecompToParent(redecomp_response_, r_gridfn);
}

void MfemMeshData::SetParentVelocity(const mfem::ParGridFunction& velocity)
{
  if (velocity_)
  {
    velocity_->SetParentGridFn(velocity);
  }
  else
  {
    velocity_ = std::make_unique<ParentField>(velocity);
  }
}

void MfemMeshData::SetLORFactor(integer lor_factor)
{
  if (lor_factor <= 1)
  {
    SLIC_WARNING_ROOT("lor_factor must be an integer > 1.  LOR factor not changed.");
    return;
  }
  if (coords_.GetParentGridFn().FESpace()->FEColl()->GetOrder() <= 1)
  {
    SLIC_WARNING_ROOT("lor_factor is only applicable to higher order geometry.  "
      "LOR factor not changed.");
    return;
  }
  lor_factor_ = lor_factor;
  // note: calls ParMesh's move ctor
  lor_mesh_ = std::make_unique<mfem::ParMesh>(mfem::ParMesh::MakeRefined(
    submesh_, lor_factor, mfem::BasisType::ClosedUniform
  ));
  submesh_lor_xfer_ = std::make_unique<SubmeshLORTransfer>(
    *submesh_xfer_gridfn_.ParFESpace(),
    *lor_mesh_
  );
}

MfemMeshData::UpdateData::UpdateData(
  mfem::ParSubMesh& submesh,
  mfem::ParMesh* lor_mesh,
  const mfem::ParFiniteElementSpace& parent_fes,
  mfem::ParGridFunction& submesh_gridfn,
  SubmeshLORTransfer* submesh_lor_xfer,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  integer num_verts_per_elem
)
: redecomp_mesh_ { lor_mesh ? 
    redecomp::RedecompMesh(*lor_mesh) :
    redecomp::RedecompMesh(submesh)
  },
  vector_xfer_ { parent_fes, submesh_gridfn, submesh_lor_xfer, redecomp_mesh_ }
{
  UpdateConnectivity(attributes_1, attributes_2, num_verts_per_elem);
}

void MfemMeshData::UpdateData::UpdateConnectivity(
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  integer num_verts_per_elem
)
{
  conn_1_.reserve(redecomp_mesh_.GetNE() * num_verts_per_elem);
  conn_2_.reserve(redecomp_mesh_.GetNE() * num_verts_per_elem);
  elem_map_1_.reserve(static_cast<size_t>(redecomp_mesh_.GetNE()));
  elem_map_2_.reserve(static_cast<size_t>(redecomp_mesh_.GetNE()));
  for (int e{}; e < redecomp_mesh_.GetNE(); ++e)
  {
    auto elem_attrib = redecomp_mesh_.GetAttribute(e);
    auto elem_conn = mfem::Array<int>();
    redecomp_mesh_.GetElementVertices(e, elem_conn);
    bool elem_on_1 = false;
    for (auto attribute_1 : attributes_1)
    {
      if (attribute_1 == elem_attrib)
      {
        elem_map_1_.push_back(e);
        conn_1_.resize(elem_map_1_.size(), num_verts_per_elem);
        for (int v{}; v < num_verts_per_elem; ++v)
        {
          conn_1_(elem_map_1_.size() - 1, v) = elem_conn[v];
        }
        elem_on_1 = true;
        break;
      }
    }
    if (!elem_on_1)
    {
      for (auto attribute_2 : attributes_2)
      {
        if (attribute_2 == elem_attrib)
        {
          elem_map_2_.push_back(e);
          conn_2_.resize(elem_map_2_.size(), num_verts_per_elem);
          for (int v{}; v < num_verts_per_elem; ++v)
          {
            conn_2_(elem_map_2_.size() - 1, v) = elem_conn[v];
          }
          break;
        }
      }
    }
  }
  conn_1_.shrink();
  conn_2_.shrink();
  elem_map_1_.shrink_to_fit();
  elem_map_2_.shrink_to_fit();
}

MfemMeshData::UpdateData& MfemMeshData::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const MfemMeshData::UpdateData& MfemMeshData::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

mfem::ParSubMesh MfemMeshData::CreateSubmesh(
  const mfem::ParMesh& parent_mesh,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2
)
{
  // TODO: Create PR for mfem::ParSubMesh::CreateFromBoundary taking a const
  // reference to attributes. Then we can construct submesh_ in the initializer
  // list without this function (because CreateFromBoundary will be willing to
  // take an rvalue for attributes)
  auto attributes_array 
    = arrayFromSet(mergeContainers(attributes_1, attributes_2));
  return mfem::ParSubMesh::CreateFromBoundary(
    parent_mesh,
    attributes_array
  );
}

MfemSubmeshData::MfemSubmeshData(
  mfem::ParSubMesh& submesh,
  mfem::ParMesh* lor_mesh,
  std::unique_ptr<mfem::FiniteElementCollection> pressure_fec,
  integer pressure_vdim
)
: submesh_pressure_ {
    new mfem::ParFiniteElementSpace(
      &submesh,
      pressure_fec.get(),
      pressure_vdim
    )
  },
  pressure_ { submesh_pressure_ },
  submesh_lor_xfer_ { lor_mesh ?
    std::make_unique<SubmeshLORTransfer>(
      *submesh_pressure_.ParFESpace(), 
      *lor_mesh
    ) :
    nullptr
  }
{
  submesh_pressure_.MakeOwner(pressure_fec.release());
  submesh_pressure_ = 0.0;
}

void MfemSubmeshData::UpdateMfemSubmeshData(redecomp::RedecompMesh& redecomp_mesh)
{
  update_data_ = std::make_unique<UpdateData>(
    *submesh_pressure_.ParFESpace(),
    submesh_lor_xfer_.get(),
    redecomp_mesh
  );
  pressure_.UpdateField(update_data_->pressure_xfer_);
  redecomp_gap_.SetSpace(pressure_.GetRedecompGridFn().FESpace());
  redecomp_gap_ = 0.0;
}

void MfemSubmeshData::GetSubmeshGap(mfem::Vector& g) const
{
  g.SetSize(submesh_pressure_.ParFESpace()->GetVSize());
  mfem::ParGridFunction g_gridfn(submesh_pressure_.ParFESpace(), g);
  g_gridfn = 0.0;
  GetPressureTransfer().RedecompToSubmesh(redecomp_gap_, g_gridfn);
}

MfemSubmeshData::UpdateData::UpdateData(
  mfem::ParFiniteElementSpace& submesh_fes,
  SubmeshLORTransfer* submesh_lor_xfer,
  redecomp::RedecompMesh& redecomp_mesh
)
: pressure_xfer_ { submesh_fes, submesh_lor_xfer, redecomp_mesh }
{}

MfemSubmeshData::UpdateData& MfemSubmeshData::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const MfemSubmeshData::UpdateData& MfemSubmeshData::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

MfemJacobianData::MfemJacobianData(
  const MfemMeshData& parent_data,
  const MfemSubmeshData& submesh_data
)
: parent_data_ { parent_data },
  submesh_data_ { submesh_data },
  block_offsets_ { 3 }
{
  SLIC_ERROR_ROOT_IF(
    parent_data.GetParentCoords().ParFESpace()->FEColl()->GetOrder() > 1,
    "Higher order meshes not yet supported for Jacobian matrices."
  );

  // NOTE: Looks like GetFrom() should be const in MFEM
  mfem::SubMeshUtils::BuildVdofToVdofMap(
    parent_data_.GetSubmeshFESpace(),
    *parent_data_.GetParentCoords().FESpace(),
    const_cast<mfem::ParSubMesh&>(parent_data_.GetSubmesh()).GetFrom(),
    parent_data_.GetSubmesh().GetParentElementIDMap(),
    submesh2parent_vdof_list_
  );

  int my_rank {0};
  MPI_Comm_rank(TRIBOL_COMM_WORLD, &my_rank);
  auto offset_idx = HYPRE_AssumedPartitionCheck() ? 0 : my_rank;
  auto dof_offset = parent_data_.GetParentCoords().ParFESpace()->
    GetDofOffsets()[offset_idx];
  for (auto& vdof : submesh2parent_vdof_list_)
  {
    vdof = vdof + dof_offset;
  }

  auto disp_size = parent_data_.GetParentCoords().ParFESpace()
    ->GetTrueVSize();
  auto lm_size = submesh_data_.GetSubmeshPressure().ParFESpace()
    ->GetTrueVSize();
  block_offsets_[0] = 0;
  block_offsets_[1] = disp_size;
  block_offsets_[2] = disp_size + lm_size;
}

void MfemJacobianData::UpdateJacobianXfer()
{
  update_data_ = std::make_unique<UpdateData>(parent_data_, submesh_data_);
}

std::unique_ptr<mfem::BlockOperator> MfemJacobianData::GetMfemBlockJacobian(
  const MethodData& method_data
) const
{
  // 0 = displacement DOFs, 1 = lagrange multiplier DOFs
  // (0,0) block is empty (for now)
  // (1,1) block is empty
  // (0,1) and (1,0) are symmetric (for now)
  const auto& elem_map_1 = parent_data_.GetElemMap1();
  const auto& elem_map_2 = parent_data_.GetElemMap2();
  auto mortar_elems = method_data
    .getBlockJElementIds()[static_cast<int>(BlockSpace::MORTAR)];
  for (auto& mortar_elem : mortar_elems)
  {
    mortar_elem = elem_map_1[static_cast<size_t>(mortar_elem)];
  }
  auto nonmortar_elems = method_data
    .getBlockJElementIds()[static_cast<int>(BlockSpace::NONMORTAR)];
  for (auto& nonmortar_elem : nonmortar_elems)
  {
    nonmortar_elem = elem_map_2[static_cast<size_t>(nonmortar_elem)];
  }
  auto lm_elems = method_data
    .getBlockJElementIds()[static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER)];
  for (auto& lm_elem : lm_elems)
  {
    lm_elem = elem_map_2[static_cast<size_t>(lm_elem)];
  }
  // get (1,0) block
  const auto& elem_J_1 = method_data.getBlockJ()(
    static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER),
    static_cast<int>(BlockSpace::MORTAR)
  );
  const auto& elem_J_2 = method_data.getBlockJ()(
    static_cast<int>(BlockSpace::LAGRANGE_MULTIPLIER),
    static_cast<int>(BlockSpace::NONMORTAR)
  );
  // move to submesh level
  auto submesh_J = GetUpdateData().submesh_redecomp_xfer_->TransferToParallelSparse(
    lm_elems, 
    mortar_elems, 
    elem_J_1
  );
  submesh_J += GetUpdateData().submesh_redecomp_xfer_->TransferToParallelSparse(
    lm_elems, 
    nonmortar_elems, 
    elem_J_2
  );
  submesh_J.Finalize();

  // transform J values from submesh to parent
  auto J = submesh_J.GetJ();
  auto submesh_vector_fes = parent_data_.GetSubmeshFESpace();
  auto mpi = redecomp::MPIUtility(submesh_vector_fes.GetComm());
  auto submesh_dof_offsets = axom::Array<int>(mpi.NRanks() + 1, mpi.NRanks() + 1);
  if (HYPRE_AssumedPartitionCheck())
  {
    submesh_dof_offsets[mpi.MyRank()+1] = submesh_vector_fes.GetDofOffsets()[1];
    mpi.Allreduce(&submesh_dof_offsets, MPI_SUM);
  }
  else
  {
    for (int i{0}; i < mpi.NRanks(); ++i)
    {
        submesh_dof_offsets[i] = submesh_vector_fes.GetDofOffsets()[i];
    }
    
  }
  auto send_J_by_rank = redecomp::MPIArray<int>(&mpi);
  auto J_idx = redecomp::MPIArray<int>(&mpi);
  auto est_num_J = submesh_J.NumNonZeroElems() / mpi.NRanks();
  for (int r{}; r < mpi.NRanks(); ++r)
  {
    if (r == mpi.MyRank())
    {
        send_J_by_rank[r].shrink();
        J_idx[r].shrink();
    }
    else
    {
        send_J_by_rank[r].reserve(est_num_J);
        J_idx[r].reserve(est_num_J);
    }
    
  }
  for (int j{}; j < submesh_J.NumNonZeroElems(); ++j)
  {
    if (J[j] >= submesh_dof_offsets[mpi.MyRank()] 
        && J[j] < submesh_dof_offsets[mpi.MyRank() + 1])
    {
        J[j] = submesh2parent_vdof_list_[J[j] - submesh_dof_offsets[mpi.MyRank()]];
    }
    else
    {
        for (int r{}; r < mpi.NRanks(); ++r)
        {
        if (J[j] >= submesh_dof_offsets[r] && J[j] < submesh_dof_offsets[r + 1])
        {
          send_J_by_rank[r].push_back(J[j] - submesh_dof_offsets[r]);
          J_idx[r].push_back(j);
          break;
        }
        }
    }
  }
  auto recv_J_by_rank = redecomp::MPIArray<int>(&mpi);
  recv_J_by_rank.SendRecvArrayEach(send_J_by_rank);
  // update recv_J_by_rank
  for (int r{}; r < mpi.NRanks(); ++r)
  {
    for (auto& recv_J : recv_J_by_rank[r])
    {
        recv_J = submesh2parent_vdof_list_[recv_J];
    }
  }
  // give it back to send_J_by_rank
  send_J_by_rank.SendRecvArrayEach(recv_J_by_rank);
  for (int r{}; r < mpi.NRanks(); ++r)
  {
    for (int j{}; j < send_J_by_rank[r].size(); ++j)
    {
        J[J_idx[r][j]] = send_J_by_rank[r][j];
    }
  }

  // create block operator
  auto block_J = std::make_unique<mfem::BlockOperator>(block_offsets_);
  block_J->owns_blocks = 1;

  // fill block operator
  auto mpi_comm = parameters_t::getInstance().problem_comm;
  auto parent_test_fes = submesh_data_.GetSubmeshPressure().ParFESpace();
  auto parent_trial_fes = parent_data_.GetParentCoords().ParFESpace();
  auto J_full = std::make_unique<mfem::HypreParMatrix>(
    mpi_comm, parent_test_fes->GetVSize(), 
    parent_test_fes->GlobalVSize(), parent_trial_fes->GlobalVSize(),
    submesh_J.GetI(), submesh_J.GetJ(), submesh_J.GetData(),
    parent_test_fes->GetDofOffsets(), parent_trial_fes->GetDofOffsets()
  );
  auto J_true = std::unique_ptr<mfem::HypreParMatrix>(mfem::RAP(
    parent_test_fes->Dof_TrueDof_Matrix(),
    J_full.get(),
    parent_trial_fes->Dof_TrueDof_Matrix()
  ));
  block_J->SetBlock(0, 1, new mfem::TransposeOperator(J_true.get()));
  block_J->SetBlock(1, 0, J_true.release());

  return block_J;
}

MfemJacobianData::UpdateData::UpdateData(
  const MfemMeshData& parent_data,
  const MfemSubmeshData& submesh_data
)
{
  auto dual_submesh_fes = &submesh_data.GetSubmeshFESpace();
  auto primal_submesh_fes = &parent_data.GetSubmeshFESpace();
  if (parent_data.GetLORMesh())
  {
    dual_submesh_fes = submesh_data.GetLORMeshFESpace();
    primal_submesh_fes = parent_data.GetLORMeshFESpace();
  }
  submesh_redecomp_xfer_ = std::make_unique<redecomp::MatrixTransfer>(
    *dual_submesh_fes,
    *primal_submesh_fes,
    *submesh_data.GetRedecompGap().FESpace(),
    *parent_data.GetRedecompResponse().FESpace()
  );
}

MfemJacobianData::UpdateData& MfemJacobianData::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const MfemJacobianData::UpdateData& MfemJacobianData::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

} // end tribol namespace

#endif /* BUILD_REDECOMP */