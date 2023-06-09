// wrapTestMesh.h
// This file is generated by Shroud 0.12.1. Do not edit.

// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

/**
 * \file wrapTestMesh.h
 * \brief Shroud generated wrapper for TestMesh class
 */
// For C users and C++ implementation

#ifndef WRAPTESTMESH_H
#define WRAPTESTMESH_H

#include "typesTRIBOL_TEST_MESH.h"

// splicer begin class.TestMesh.CXX_declarations
// splicer end class.TestMesh.CXX_declarations

#ifdef __cplusplus
extern "C" {
#endif

// splicer begin class.TestMesh.C_declarations
// splicer end class.TestMesh.C_declarations

TRIBOL_TEST_MESH_TestMesh * TRIBOL_TEST_MESH_TestMesh_new(TRIBOL_TEST_MESH_TestMesh * SHC_rv);

void TRIBOL_TEST_MESH_TestMesh_delete(TRIBOL_TEST_MESH_TestMesh * self);

void TRIBOL_TEST_MESH_TestMesh_setup_contact_mesh_hex(TRIBOL_TEST_MESH_TestMesh * self, int numElemsX1, int numElemsY1, int numElemsZ1, double xMin1, double yMin1, double zMin1, double xMax1, double yMax1, double zMax1, int numElemsX2, int numElemsY2, int numElemsZ2, double xMin2, double yMin2, double zMin2, double xMax2, double yMax2, double zMax2, double thetaMortar, double thetaNonmortar);

double * TRIBOL_TEST_MESH_TestMesh_get_x(const TRIBOL_TEST_MESH_TestMesh * self);

double * TRIBOL_TEST_MESH_TestMesh_get_x_bufferify(const TRIBOL_TEST_MESH_TestMesh * self, TRIBOL_TEST_MESH_SHROUD_array *DSHC_rv);

double * TRIBOL_TEST_MESH_TestMesh_get_y(const TRIBOL_TEST_MESH_TestMesh * self);

double * TRIBOL_TEST_MESH_TestMesh_get_y_bufferify(const TRIBOL_TEST_MESH_TestMesh * self, TRIBOL_TEST_MESH_SHROUD_array *DSHC_rv);

double * TRIBOL_TEST_MESH_TestMesh_get_z(const TRIBOL_TEST_MESH_TestMesh * self);

double * TRIBOL_TEST_MESH_TestMesh_get_z_bufferify(const TRIBOL_TEST_MESH_TestMesh * self, TRIBOL_TEST_MESH_SHROUD_array *DSHC_rv);

int TRIBOL_TEST_MESH_TestMesh_get_mortar_face_connectivity_size(const TRIBOL_TEST_MESH_TestMesh * self);

int * TRIBOL_TEST_MESH_TestMesh_get_mortar_face_connectivity(const TRIBOL_TEST_MESH_TestMesh * self);

int * TRIBOL_TEST_MESH_TestMesh_get_mortar_face_connectivity_bufferify(const TRIBOL_TEST_MESH_TestMesh * self, TRIBOL_TEST_MESH_SHROUD_array *DSHC_rv);

int TRIBOL_TEST_MESH_TestMesh_get_nonmortar_face_connectivity_size(const TRIBOL_TEST_MESH_TestMesh * self);

int * TRIBOL_TEST_MESH_TestMesh_get_nonmortar_face_connectivity(const TRIBOL_TEST_MESH_TestMesh * self);

int * TRIBOL_TEST_MESH_TestMesh_get_nonmortar_face_connectivity_bufferify(const TRIBOL_TEST_MESH_TestMesh * self, TRIBOL_TEST_MESH_SHROUD_array *DSHC_rv);

int TRIBOL_TEST_MESH_TestMesh_get_numtotalnodes(TRIBOL_TEST_MESH_TestMesh * self);

int TRIBOL_TEST_MESH_TestMesh_get_nummortarfaces(TRIBOL_TEST_MESH_TestMesh * self);

int TRIBOL_TEST_MESH_TestMesh_get_numnonmortarfaces(TRIBOL_TEST_MESH_TestMesh * self);

#ifdef __cplusplus
}
#endif

#endif  // WRAPTESTMESH_H
