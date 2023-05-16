// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "tribol/mesh/MfemData.hpp"

#include "axom/slic.hpp"

namespace tribol
{

SubmeshRedecompTransfer::SubmeshRedecompTransfer(
  redecomp::RedecompMesh& redecomp,
  const mfem::ParFiniteElementSpace& submesh_fes
)
: submesh_fes_ { submesh_fes },
  redecomp_fes_ { std::make_unique<mfem::FiniteElementSpace>(
    &redecomp,
    submesh_fes_.FEColl(),
    redecomp.SpaceDimension(),
    mfem::Ordering::byNODES
  ) },
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
  auto dst = mfem::GridFunction(redecomp_fes_.get());
  redecomp_xfer_.TransferToSerial(src, dst);
  return dst;
}

mfem::ParGridFunction SubmeshRedecompTransfer::RedecompToSubmesh(
  const mfem::GridFunction& src
) const
{
  auto dst = mfem::ParGridFunction(&submesh_fes_);
  redecomp_xfer_.TransferToParallel(src, dst);
  return dst;
}

void SubmeshRedecompTransfer::RedecompToSubmesh(
  const mfem::GridFunction& src,
  mfem::ParGridFunction& dst
) const
{
  redecomp_xfer_.TransferToParallel(src, dst);
}

void SubmeshRedecompTransfer::UpdateRedecomp(redecomp::RedecompMesh& redecomp)
{
  redecomp_fes_ = std::make_unique<mfem::FiniteElementSpace>(
    &redecomp,
    submesh_fes_.FEColl(),
    redecomp.SpaceDimension(),
    mfem::Ordering::byNODES
  );
}

ParentRedecompTransfer::ParentRedecompTransfer(
  redecomp::RedecompMesh& redecomp,
  mfem::ParSubMesh& submesh,
  const mfem::ParFiniteElementSpace& parent_fes
)
: parent_fes_ { parent_fes },
  redecomp_xfer_ { redecomp, mfem::ParFiniteElementSpace(
    &submesh,
    parent_fes_.FEColl(),
    submesh.SpaceDimension(),
    mfem::Ordering::byNODES
  ) },
  submesh_gridfn_ { redecomp_xfer_.EmptySubmeshGridFn() }
{
  SLIC_ERROR_ROOT_IF(
    submesh.GetParent() != parent_fes_.GetParMesh(),
    "submesh's parent must match the parent_fes ParMesh."
  );
}

mfem::GridFunction ParentRedecompTransfer::ParentToRedecomp(
  const mfem::ParGridFunction& src
) const
{
  redecomp_xfer_.GetSubmesh().Transfer(src, submesh_gridfn_);
  return redecomp_xfer_.SubmeshToRedecomp(submesh_gridfn_);
}

mfem::ParGridFunction ParentRedecompTransfer::RedecompToParent(
  const mfem::GridFunction& src
) const
{
  auto dst = mfem::ParGridFunction(&parent_fes_);
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

std::vector<const real*> PrimalField::GetFieldPtrs() const
{
  auto data_ptrs = std::vector<const real*>(3, nullptr);
  for (size_t i{}; i < static_cast<size_t>(GetRedecompField().FESpace()->GetVDim()); ++i)
  {
    data_ptrs[i] = &GetRedecompField()(GetRedecompField().FESpace()->DofToVDof(0, i));
  }
  return data_ptrs;
}

std::vector<real*> PrimalField::GetFieldPtrs(mfem::GridFunction& redecomp_fn)
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

MfemMeshData::MfemMeshData(
  integer mesh_id_1,
  integer mesh_id_2,
  mfem::ParMesh& mesh,
  const mfem::ParGridFunction& current_coords,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  std::unique_ptr<const mfem::FiniteElementCollection> dual_fec,
  integer dual_vdim
)
: mesh_id_1_ { mesh_id_1 },
  mesh_id_2_ { mesh_id_2 },
  mesh_ { mesh },
  attributes_1_ { attributes_1 },
  attributes_2_ { attributes_2 },
  submesh_ { CreateSubmesh(mesh_, attributes_1_, attributes_2_) },
  coords_ { current_coords },
  dual_fec_ { std::move(dual_fec) },
  dual_vdim_ { dual_vdim }
{
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
    attributes_1_, 
    attributes_2_, 
    num_verts_per_elem_,
    *coords_.GetParentField().ParFESpace(),
    dual_fec_.get(),
    dual_vdim_
  );
  coords_.UpdateField(update_data_->primal_xfer_);
  response_ = mfem::GridFunction(coords_.GetRedecompField().FESpace());
}

MfemMeshData::UpdateData::UpdateData(
  mfem::ParSubMesh& submesh,
  const std::set<integer>& attributes_1,
  const std::set<integer>& attributes_2,
  integer num_verts_per_elem,
  const mfem::ParFiniteElementSpace& parent_fes,
  const mfem::FiniteElementCollection* dual_fec,
  integer dual_vdim
)
: redecomp_ { submesh },
  primal_xfer_ { redecomp_, submesh, parent_fes },
  dual_xfer_ { dual_fec != nullptr ?
    std::make_unique<SubmeshRedecompTransfer>(redecomp_, mfem::ParFiniteElementSpace(
      &submesh,
      dual_fec,
      dual_vdim,
      mfem::Ordering::byNODES
    )) :
    nullptr
  }
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

} // end tribol namespace
