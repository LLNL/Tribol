// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/MfemData.hpp"

#include "axom/slic.hpp"

namespace tribol
{

SubmeshLORTransfer::SubmeshLORTransfer(
  mfem::ParFiniteElementSpace& ho_submesh_fes,
  mfem::ParMesh& lor_submesh
)
: lor_gridfn_ { CreateLORGridFunction(
    lor_submesh,
    std::make_unique<mfem::H1_FECollection>(1, lor_submesh.SpaceDimension()),
    ho_submesh_fes.GetVDim()
  ) },
  lor_xfer_ { ho_submesh_fes, *lor_gridfn_.ParFESpace() }
{}

void SubmeshLORTransfer::TransferToLORGridFn(
  const mfem::ParGridFunction& ho_src
)
{
  lor_xfer_.ForwardOperator().Mult(ho_src, lor_gridfn_);
}

void SubmeshLORTransfer::TransferFromLORGridFn(
  mfem::ParGridFunction& ho_dst
) const
{
  lor_xfer_.BackwardOperator().Mult(lor_gridfn_, ho_dst);
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
  SubmeshLORTransfer* lor_xfer,
  redecomp::RedecompMesh& redecomp
)
: submesh_fes_ { submesh_fes },
  redecomp_fes_ { CreateRedecompFESpace(redecomp, submesh_fes_) },
  lor_xfer_ { lor_xfer },
  redecomp_xfer_ { } // default (element transfer) constructor
{
  // make sure submesh_fes is a submesh and redecomp's parent is submesh_fes's
  // submesh
  SLIC_ERROR_ROOT_IF(
    !mfem::ParSubMesh::IsParSubMesh(submesh_fes_.GetParMesh()),
    "submesh_fes must be on a ParSubMesh."
  );
  SLIC_ERROR_ROOT_IF(
    &redecomp.getParent() != submesh_fes_.GetParMesh(),
    "redecomp's parent must match the submesh_fes ParMesh."
  );
}

mfem::GridFunction SubmeshRedecompTransfer::SubmeshToRedecomp(
  const mfem::ParGridFunction& src
) const
{
  auto src_ptr = &src;
  auto dst = mfem::GridFunction(redecomp_fes_.get());
  dst = 0.0;
  if (lor_xfer_)
  {
    lor_xfer_->GetLORGridFn() = 0.0;
    lor_xfer_->TransferToLORGridFn(src);
    src_ptr = &lor_xfer_->GetLORGridFn();
  }
  redecomp_xfer_.TransferToSerial(*src_ptr, dst);
  return dst;
}

mfem::ParGridFunction SubmeshRedecompTransfer::RedecompToSubmesh(
  const mfem::GridFunction& src
) const
{
  auto dst = mfem::ParGridFunction(&submesh_fes_);
  dst = 0.0;
  RedecompToSubmesh(src, dst);
  return dst;
}

void SubmeshRedecompTransfer::RedecompToSubmesh(
  const mfem::GridFunction& src,
  mfem::ParGridFunction& dst
) const
{
  auto dst_ptr = &dst;
  if (lor_xfer_)
  {
    lor_xfer_->GetLORGridFn() = 0.0;
    dst_ptr = &lor_xfer_->GetLORGridFn();
  }
  redecomp_xfer_.TransferToParallel(src, *dst_ptr);
  if (lor_xfer_)
  {
    lor_xfer_->TransferFromLORGridFn(dst);
  }
}

void SubmeshRedecompTransfer::UpdateRedecomp(redecomp::RedecompMesh& redecomp)
{
  redecomp_fes_ = CreateRedecompFESpace(redecomp, submesh_fes_);
}

std::unique_ptr<mfem::FiniteElementSpace> SubmeshRedecompTransfer::CreateRedecompFESpace(
  redecomp::RedecompMesh& redecomp,
  mfem::ParFiniteElementSpace& submesh_fes
)
{
  return std::make_unique<mfem::FiniteElementSpace>(
    &redecomp,
    submesh_fes.FEColl(),
    submesh_fes.GetVDim(),
    mfem::Ordering::byNODES
  );
}

ParentRedecompTransfer::ParentRedecompTransfer(
  const mfem::ParFiniteElementSpace& parent_fes,
  mfem::ParGridFunction& submesh_gridfn,
  SubmeshLORTransfer* lor_xfer,
  redecomp::RedecompMesh& redecomp
)
: parent_fes_ { parent_fes },
  submesh_gridfn_ { submesh_gridfn },
  redecomp_xfer_ { 
    *submesh_gridfn_.ParFESpace(),
    lor_xfer, 
    redecomp
  }
{
  // Note: this is checked in the SubmeshRedecompTransfer constructor
  // SLIC_ERROR_ROOT_IF(
  //   !mfem::ParSubMesh::IsParSubMesh(submesh_gridfn_.ParFESpace()->GetParMesh()),
  //   "submesh_gridfn_ must be associated with an mfem::ParSubMesh."
  // );
  SLIC_ERROR_ROOT_IF(
    redecomp_xfer_.GetSubmesh().GetParent() != parent_fes_.GetParMesh(),
    "submesh_gridfn's parent mesh must match the parent_fes ParMesh."
  );
}

mfem::GridFunction ParentRedecompTransfer::ParentToRedecomp(
  const mfem::ParGridFunction& src
) const
{
  submesh_gridfn_ = 0.0;
  redecomp_xfer_.GetSubmesh().Transfer(src, submesh_gridfn_);
  return redecomp_xfer_.SubmeshToRedecomp(submesh_gridfn_);
}

mfem::ParGridFunction ParentRedecompTransfer::RedecompToParent(
  const mfem::GridFunction& src
) const
{
  auto dst = mfem::ParGridFunction(&parent_fes_);
  dst = 0.0;
  submesh_gridfn_ = 0.0;
  redecomp_xfer_.RedecompToSubmesh(src, submesh_gridfn_);
  redecomp_xfer_.GetSubmesh().Transfer(submesh_gridfn_, dst);
  return dst;
}

void ParentRedecompTransfer::UpdateRedecomp(redecomp::RedecompMesh& redecomp)
{
  redecomp_xfer_.UpdateRedecomp(redecomp);
}

PrimalField::PrimalField(
  const mfem::ParGridFunction& parent
)
: parent_ { parent }
{}

void PrimalField::SetParentField(const mfem::ParGridFunction& parent)
{
  parent_ = parent;
}

void PrimalField::UpdateField(const ParentRedecompTransfer& xfer)
{
  update_data_ = std::make_unique<UpdateData>(xfer, parent_);
}

std::vector<const real*> PrimalField::GetRedecompFieldPtrs() const
{
  auto data_ptrs = std::vector<const real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(GetRedecompField().FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &GetRedecompField()(GetRedecompField().FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

std::vector<real*> PrimalField::GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn)
{
  auto data_ptrs = std::vector<real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(redecomp_fn.FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &redecomp_fn(redecomp_fn.FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

PrimalField::UpdateData& PrimalField::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const PrimalField::UpdateData& PrimalField::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

PrimalField::UpdateData::UpdateData(
  const ParentRedecompTransfer& xfer,
  const mfem::ParGridFunction& parent
)
: xfer_ { xfer },
  redecomp_ { xfer_.ParentToRedecomp(parent) }
{}

PressureField::PressureField(
  const mfem::ParGridFunction& submesh
)
: submesh_ { submesh }
{}

void PressureField::SetSubmeshField(const mfem::ParGridFunction& submesh)
{
  submesh_ = submesh;
}

void PressureField::UpdateField(const SubmeshRedecompTransfer& xfer)
{
  update_data_ = std::make_unique<UpdateData>(xfer, submesh_);
}

std::vector<const real*> PressureField::GetRedecompFieldPtrs() const
{
  auto data_ptrs = std::vector<const real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(GetRedecompField().FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &GetRedecompField()(GetRedecompField().FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

std::vector<real*> PressureField::GetRedecompFieldPtrs(mfem::GridFunction& redecomp_fn)
{
  auto data_ptrs = std::vector<real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(redecomp_fn.FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &redecomp_fn(redecomp_fn.FESpace()->DofToVDof(0, i));
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
  const SubmeshRedecompTransfer& xfer,
  const mfem::ParGridFunction& submesh
)
: xfer_ { xfer },
  redecomp_ { xfer_.SubmeshToRedecomp(submesh) }
{}

MfemMeshData::MfemMeshData(
  integer mesh_id_1,
  integer mesh_id_2,
  mfem::ParMesh& mesh,
  const mfem::ParGridFunction& current_coords,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2
)
: mesh_id_1_ { mesh_id_1 },
  mesh_id_2_ { mesh_id_2 },
  mesh_ { mesh },
  attributes_1_ { attributes_1 },
  attributes_2_ { attributes_2 },
  submesh_ { CreateSubmesh(mesh_, attributes_1_, attributes_2_) },
  coords_ { current_coords },
  submesh_xfer_gridfn_ { 
    CreateSubmeshGridFn(submesh_, *current_coords.ParFESpace())
  },
  lor_factor_ { 0 }
{
  // build LOR submesh
  if (current_coords.FESpace()->FEColl()->GetOrder() > 1)
  {
    SetLowOrderRefinedFactor(current_coords.FESpace()->FEColl()->GetOrder());
  }

  // set the element type
  mfem::Element::Type element_type = mfem::Element::QUADRILATERAL;
  if (submesh_.GetNE() > 0)
  {
    element_type = submesh_.GetElementType(0);
  }

  switch (element_type) 
  {
    case mfem::Element::POINT:
      elem_type_ = NODE;
      break;
    case mfem::Element::SEGMENT:
      elem_type_ = EDGE;
      break;
    case mfem::Element::QUADRILATERAL:
      elem_type_ = FACE;
      break;
    case mfem::Element::HEXAHEDRON:
      elem_type_ = CELL;
      break;

    case mfem::Element::TRIANGLE:
    case mfem::Element::TETRAHEDRON:
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
  coords_.SetParentField(current_coords);
}

void MfemMeshData::UpdateMeshData()
{
  update_data_ = std::make_unique<UpdateData>(
    submesh_,
    lor_submesh_.get(),
    *coords_.GetParentField().ParFESpace(),
    submesh_xfer_gridfn_,
    lor_xfer_.get(),
    attributes_1_, 
    attributes_2_, 
    num_verts_per_elem_
  );
  coords_.UpdateField(update_data_->primal_xfer_);
  redecomp_response_.SetSpace(coords_.GetRedecompField().FESpace());
  redecomp_response_ = 0.0;
  if (velocity_)
  {
    velocity_->UpdateField(update_data_->primal_xfer_);
  }
}

void MfemMeshData::SetParentVelocity(const mfem::ParGridFunction& velocity)
{
  if (velocity_)
  {
    velocity_->SetParentField(velocity);
  }
  else
  {
    velocity_ = std::make_unique<PrimalField>(velocity);
  }
}

void MfemMeshData::SetLowOrderRefinedFactor(integer lor_factor)
{
  if (lor_factor <= 1)
  {
    SLIC_WARNING_ROOT("lor_factor must be an integer > 1.  LOR factor not changed.");
    return;
  }
  if (coords_.GetParentField().FESpace()->FEColl()->GetOrder() <= 1)
  {
    SLIC_WARNING_ROOT("lor_factor is only applicable to higher order geometry.  "
      "LOR factor not changed.");
    return;
  }
  lor_factor_ = lor_factor;
  // note: calls ParMesh's move ctor
  lor_submesh_ = std::make_unique<mfem::ParMesh>(mfem::ParMesh::MakeRefined(
    submesh_, lor_factor, mfem::BasisType::ClosedUniform
  ));
  lor_xfer_ = std::make_unique<SubmeshLORTransfer>(
    *submesh_xfer_gridfn_.ParFESpace(),
    *lor_submesh_
  );
}

MfemMeshData::UpdateData::UpdateData(
  mfem::ParSubMesh& submesh,
  mfem::ParMesh* lor_submesh,
  const mfem::ParFiniteElementSpace& parent_fes,
  mfem::ParGridFunction& submesh_gridfn,
  SubmeshLORTransfer* lor_xfer,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  integer num_verts_per_elem
)
: redecomp_ { lor_submesh ? 
    redecomp::RedecompMesh(submesh) :
    redecomp::RedecompMesh(*lor_submesh)
  },
  primal_xfer_ { parent_fes, submesh_gridfn, lor_xfer, redecomp_ }
{
  UpdateConnectivity(attributes_1, attributes_2, num_verts_per_elem);
}

void MfemMeshData::UpdateData::UpdateConnectivity(
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  integer num_verts_per_elem
)
{
  conn_1_.reserve(redecomp_.GetNE() * num_verts_per_elem);
  conn_2_.reserve(redecomp_.GetNE() * num_verts_per_elem);
  elem_map_1_.reserve(static_cast<size_t>(redecomp_.GetNE()));
  elem_map_2_.reserve(static_cast<size_t>(redecomp_.GetNE()));
  for (int e{}; e < redecomp_.GetNE(); ++e)
  {
    auto elem_attrib = redecomp_.GetAttribute(e);
    auto elem_conn = mfem::Array<int>();
    redecomp_.GetElementVertices(e, elem_conn);
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
  const mfem::ParMesh& parent,
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
    parent,
    attributes_array
  );
}

mfem::ParGridFunction MfemMeshData::CreateSubmeshGridFn(
  mfem::ParSubMesh& submesh,
  mfem::ParFiniteElementSpace& parent_fes
)
{
  std::unique_ptr<mfem::FiniteElementCollection> submesh_fec {
    parent_fes.FEColl()->Clone(parent_fes.FEColl()->GetOrder())
  };
  mfem::ParGridFunction submesh_gridfn {
    new mfem::ParFiniteElementSpace(
      &submesh,
      submesh_fec.get(),
      parent_fes.GetVDim(),
      mfem::Ordering::byNODES
    )
  };
  submesh_gridfn.MakeOwner(submesh_fec.release());
  return submesh_gridfn;
}

MfemDualData::MfemDualData(
  mfem::ParSubMesh& submesh,
  mfem::ParMesh* lor_submesh,
  std::unique_ptr<mfem::FiniteElementCollection> dual_fec,
  integer dual_vdim
)
: submesh_pressure_ {
    new mfem::ParFiniteElementSpace(
      &submesh,
      dual_fec.get(),
      dual_vdim
    )
  },
  pressure_ { submesh_pressure_ },
  lor_xfer_ { lor_submesh ?
    std::make_unique<SubmeshLORTransfer>(
      *submesh_pressure_.ParFESpace(), 
      *lor_submesh
    ) :
    nullptr
  }
{
  submesh_pressure_.MakeOwner(dual_fec.release());
  submesh_pressure_ = 0.0;
}

void MfemDualData::UpdateDualData(redecomp::RedecompMesh& redecomp)
{
  update_data_ = std::make_unique<UpdateData>(
    *submesh_pressure_.ParFESpace(),
    lor_xfer_.get(),
    redecomp
  );
  pressure_.UpdateField(update_data_->dual_xfer_);
  redecomp_gap_.SetSpace(pressure_.GetRedecompField().FESpace());
  redecomp_gap_ = 0.0;
}

MfemDualData::UpdateData::UpdateData(
  mfem::ParFiniteElementSpace& submesh_fes,
  SubmeshLORTransfer* lor_xfer,
  redecomp::RedecompMesh& redecomp
)
: dual_xfer_ { submesh_fes, lor_xfer, redecomp }
{}

MfemDualData::UpdateData& MfemDualData::GetUpdateData()
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

const MfemDualData::UpdateData& MfemDualData::GetUpdateData() const
{
  SLIC_ERROR_ROOT_IF(
    update_data_ == nullptr,
    "UpdateField() must be called to generate UpdateData."
  );
  return *update_data_;
}

} // end tribol namespace
