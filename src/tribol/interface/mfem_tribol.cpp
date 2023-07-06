// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#include "mfem_tribol.hpp"

#ifdef BUILD_REDECOMP

// Tribol includes
#include "tribol.hpp"
#include "tribol/mesh/CouplingSchemeManager.hpp"

namespace tribol
{

void registerMfemMesh( integer cs_id,
                      integer mesh_id_1,
                      integer mesh_id_2,
                      mfem::ParMesh& mesh,
                      const mfem::ParGridFunction& current_coords,
                      const std::set<integer>& attributes_1,
                      const std::set<integer>& attributes_2,
                      integer contact_mode,
                      integer contact_case,
                      integer contact_method,
                      integer contact_model,
                      integer enforcement_method,
                      integer binning_method)
{
   std::unique_ptr<mfem::FiniteElementCollection> dual_fec = nullptr;
   integer dual_vdim = 0;
   auto mfem_data = std::make_unique<MfemMeshData>(
      mesh_id_1,
      mesh_id_2,
      mesh,
      current_coords,
      attributes_1,
      attributes_2
   );
   // register empty meshes so the coupling scheme is valid
   registerMesh(
      mesh_id_1, 0, 0, nullptr, mfem_data->GetElemType(), nullptr, nullptr);
   registerMesh(
      mesh_id_2, 0, 0, nullptr, mfem_data->GetElemType(), nullptr, nullptr);
   registerCouplingScheme(
      cs_id,
      mesh_id_1,
      mesh_id_2,
      contact_mode,
      contact_case,
      contact_method,
      contact_model,
      enforcement_method,
      binning_method
   );
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   if (enforcement_method == LAGRANGE_MULTIPLIER)
   {
      dual_fec = std::make_unique<mfem::H1_FECollection>(
         current_coords.FESpace()->FEColl()->GetOrder(),
         mesh.SpaceDimension()
      );
      if (contact_model == FRICTIONLESS)
      {
         dual_vdim = 1;
      }
      else if (contact_model == TIED || contact_model == COULOMB)
      {
         dual_vdim = mesh.SpaceDimension();
      }
      else
      {
         SLIC_ERROR_ROOT("Unsupported contact model. "
           "Only FRICTIONLESS, TIED, and COULOMB supported.");
      }
      coupling_scheme->setMfemSubmeshData(
         std::make_unique<MfemSubmeshData>(
            mfem_data->GetSubmesh(),
            mfem_data->GetLORMesh(),
            std::move(dual_fec),
            dual_vdim
         )
      );
      auto lm_options = coupling_scheme->getEnforcementOptions().lm_implicit_options;
      if (
         lm_options.enforcement_option_set && 
         (
            lm_options.eval_mode == ImplicitEvalMode::MORTAR_JACOBIAN ||
            lm_options.eval_mode == ImplicitEvalMode::MORTAR_RESIDUAL_JACOBIAN
         )
      )
      {
         coupling_scheme->setMatrixXfer(std::make_unique<MfemJacobianData>(
            *mfem_data,
            *coupling_scheme->getMfemSubmeshData()
         ));
      }
   }
   coupling_scheme->setMfemMeshData(std::move(mfem_data));

}

void setMfemLowOrderRefinedFactor( integer cs_id,
                                   integer lor_factor)
{
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemData(),
      "Coupling scheme does not contain MFEM data. "
      "Was the coupling scheme created with registerParMesh()?"
   );
   coupling_scheme->getMfemMeshData()->SetLORFactor(lor_factor);
}

void registerMfemVelocity( integer cs_id, const mfem::ParGridFunction& v )
{
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemData(), 
      "Coupling scheme does not contain MFEM data. "
      "Was the coupling scheme created with registerParMesh()?"
   );
   coupling_scheme->getMfemMeshData()->SetParentVelocity(v);
}

void getMfemResponse( integer cs_id, mfem::ParGridFunction& r )
{
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemData(), 
      "Coupling scheme does not contain MFEM data. "
      "Was the coupling scheme created with registerParMesh()?"
   );
   coupling_scheme->getMfemMeshData()->GetParentResponse(r);
}

std::unique_ptr<mfem::BlockOperator> getMfemBlockJacobian( integer csId )
{
   CouplingScheme* coupling_scheme = CouplingSchemeManager::getInstance().
      getCoupling(csId);
   SparseMode sparse_mode = coupling_scheme
      ->getEnforcementOptions().lm_implicit_options.sparse_mode;
   if (sparse_mode != SparseMode::MFEM_ELEMENT_DENSE)
   {
      SLIC_ERROR_ROOT("Jacobian is assembled and can be accessed by " 
         "getMfemSparseMatrix() or getCSRMatrix(). For (unassembled) element "
         "Jacobian contributions, call setLagrangeMultiplierOptions() with "
         "SparseMode::MFEM_ELEMENT_DENSE before calling update().");
   }
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemData(), 
      "No MFEM data exists. The coupling scheme "
      "must be registered using registerParMesh() to use this method."
   );
   return coupling_scheme->getMfemJacobianData()->GetMfemBlockJacobian(
      *coupling_scheme->getMethodData()
   );
}

void getMfemGap( integer cs_id, mfem::ParGridFunction& g )
{
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemSubmeshData(), 
      "Coupling scheme does not contain MFEM dual field data. "
      "Was the coupling scheme created with registerParMesh() and is the "
      "enforcement method LAGRANGE_MULTIPLIER?"
   );
   coupling_scheme->getMfemSubmeshData()->GetSubmeshGap(g);
}

mfem::ParGridFunction& getMfemPressure( integer cs_id )
{
   auto coupling_scheme = CouplingSchemeManager::getInstance().getCoupling(cs_id);
   SLIC_ERROR_ROOT_IF(
      !coupling_scheme->hasMfemSubmeshData(), 
      "Coupling scheme does not contain MFEM dual field data. "
      "Was the coupling scheme created with registerParMesh() and is the "
      "enforcement method LAGRANGE_MULTIPLIER?"
   );
   return coupling_scheme->getMfemSubmeshData()->GetSubmeshPressure();
}

void updateMfemParallelDecomposition()
{
   CouplingSchemeManager& csManager = CouplingSchemeManager::getInstance();
   int numCouplings = csManager.getNumberOfCouplings();

   for(int csIndex =0; csIndex < numCouplings; ++csIndex)
   {
      if(!csManager.hasCoupling(csIndex))
      {
         continue;
      }

      CouplingScheme* couplingScheme  = csManager.getCoupling(csIndex);

      // update redecomp meshes if supplied mfem data
      if (couplingScheme->hasMfemData())
      {
         auto mfem_data = couplingScheme->getMfemMeshData();
         axom::Array<int> mesh_ids {2, 2};
         mesh_ids[0] = mfem_data->GetMesh1ID();
         mesh_ids[1] = mfem_data->GetMesh2ID();
         // creates a new RedecompMesh and updates grid functions on RedecompMesh
         mfem_data->UpdateMeshData();
         auto coord_ptrs = mfem_data->GetRedecompCoordsPtrs();

         registerMesh(
            mesh_ids[0],
            mfem_data->GetMesh1NE(),
            mfem_data->GetNV(),
            mfem_data->GetMesh1Conn(),
            mfem_data->GetElemType(),
            coord_ptrs[0],
            coord_ptrs[1],
            coord_ptrs[2]
         );
         registerMesh(
            mesh_ids[1],
            mfem_data->GetMesh2NE(),
            mfem_data->GetNV(),
            mfem_data->GetMesh2Conn(),
            mfem_data->GetElemType(),
            coord_ptrs[0],
            coord_ptrs[1],
            coord_ptrs[2]
         );

         auto f_ptrs = mfem_data->GetRedecompResponsePtrs();
         registerNodalResponse(
            mesh_ids[0], f_ptrs[0], f_ptrs[1], f_ptrs[2]);
         registerNodalResponse(
            mesh_ids[1], f_ptrs[0], f_ptrs[1], f_ptrs[2]);
         if (mfem_data->HasVelocity())
         {
            auto v_ptrs = mfem_data->GetRedecompVelocityPtrs();
            registerNodalVelocities(
               mesh_ids[0], v_ptrs[0], v_ptrs[1], v_ptrs[2]);
            registerNodalVelocities(
               mesh_ids[1], v_ptrs[0], v_ptrs[1], v_ptrs[2]);
         }
         if (couplingScheme->getEnforcementMethod() == LAGRANGE_MULTIPLIER)
         {
            SLIC_ERROR_ROOT_IF(couplingScheme->getContactModel() != FRICTIONLESS,
              "Only frictionless contact is supported.");
            auto dual_data = couplingScheme->getMfemSubmeshData();
            dual_data->UpdateSubmeshData(mfem_data->GetRedecompMesh());
            auto g_ptrs = dual_data->GetRedecompGapPtrs();
            registerMortarGaps(mesh_ids[1], g_ptrs[0]);
            auto p_ptrs = dual_data->GetRedecompPressurePtrs();
            registerMortarPressures(mesh_ids[1], p_ptrs[0]);
            auto lm_options = couplingScheme->getEnforcementOptions().lm_implicit_options;
            if (couplingScheme->hasMfemJacobianData())
            {
               couplingScheme->getMfemJacobianData()->UpdateJacobianXfer();
            }
         }
      }

   }

}

} // end tribol namespace

#endif /* BUILD_REDECOMP */
