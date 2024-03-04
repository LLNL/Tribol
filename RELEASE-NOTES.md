# TRIBOL -- Interface Physics Library -- Release Notes

All notable changes to this project will be documented in this file

The format of this file is based on [Keep a
Changelog](http://keepachangelog.com/en/1.0.0/).

## [Version 0.2.0] - Release date YYYY-MM-DD

### Added
- Added support for a linear tetrahedral mesh under the TestMesh class.
- Added coupling scheme tests with null pointer registration.
- Multi-rank contact API using MFEM data structures.
- New API calls for MFEM data structures (see `interface/mfem_tribol.hpp`).
- Updated the penalty stiffness calculation using the `KINEMATIC_CONSTANT` option
  to follow the `springs-in-serial` stiffness model used for `KINEMATIC_ELEMENT`.
- API function to support turning on or off the timestep calculation for 
  the common plane method.

### Changed
- Return negative timestep vote for non-null meshes with null velocity pointers.
- Refactored how surface elements are characterized, now using finite element
  order and type of cell.
- Renamed getMfemSparseMatrix() to getJacobianSparseMatrix() and getCSRMatrix()
  to getJacobianCSRMatrix() to avoid confusion with the separate new MFEM
  interface.
- Logging refactor using SLIC macros. Lots of warnings were demoted to DEBUG level.
- Changed various computational geometry routines to return a FaceGeomError enum error handling
- Changed LoggingLevel enum names by appending a TRIBOL prefix to avoid MACRO conflicts
  with host codes.
- Removed nullptr errors to allow more function call uniformity for ranks with null meshes. 
  Also removed any `continue` statements for null meshes.
- Updated logging in timestep vote by removing logging macro calls inside the interface pairs 
  loop; don't error out in the presence of a bad dt vote, but issue debug print.
  
### Fixed
- Allow null velocity and response pointers for various use cases
- Tolerancing bug that produced negative timestep estimates in the presence of numerically
  zero face velocities.

## [Version 0.1.0] - Release date 2023-04-21

### Added
- Licensing information for open source release.
- NOTICE file in root directory.
- LICENSE file in root directory.
- Binning for intermediate level search utilizing Axom.
- Common plane with single integration point penalty enforcement.
- Single mortar method per Puso 2004 with Lagrange multiplier enforcment. Note,
  Jacobian is in simplified, symmetric form.
- Mortar Jacobian output in linked list, CSR, or unassembled element formats.
- Aligned mortar method that assumes perfect, conforming alignment of each
  opposing surface face with no sliding.
- Mortar weights with CSR output.
- API function to set plot options for VTK visualization.
- Tied contact model with common plane method and penalty enforcement.
- Element-based penalty stiffness calculation can be used when arrays of element
  thicknesses and bulk moduli are registered.
- Added/updated ways to register element fields for the element-based penalty.
- Timestep vote based on gap and velocity projection.
- Added register routines for nodal and element based fields for integers and
  reals.
- Added support for gap-rate penalization for the common plane method. See
  setPenaltyOptions().

### Changed
- Examples now use `CLI11` for command line argument parsing.
- Improves tribol's build system, including tribol's import targets 
  (for importing tribol into other projects).
- Renamed the `TPL` directory to `external` and improves how internal TPLs are
  built.
- Updates `Axom` submodule to commit 2fc51ec (post version 0.3.3).
- Updates `blt` submodule to version 0.3.0.
- `m_numNodes` on the MeshData class is changed to `m_lengthNodalData`. This is
  strictly the number of array elements in the registered nodal data, such as
  nodal positions and nodal forces.
- Flipped the COMMON_PLANE normal definition to be commensurate with the gap
  convention set forth by the mortar implementation. Had to flip the summation
  (+/-) in the contact nodal forces for the COMMON_PLANE method.
- Added a data directory to hold externally generated meshes used in some new
  tests.
- Consolidated visualization routines that write to VTK for contact plane
  methods.
- No longer take ratios of bulk-modulus / element-thickness, but instead take
  arrays of those quantities separately to use in both element-based penalty
  calculations and the timestep vote. 
- setPenaltyStiffness() has been deprecated in lieu of a new function,
  setPenaltyData(), which is used to specify "single" or "element" penalty and
  to register any required data.
- The previous generic field registration routines (e.g. registerNodalField)
  have been modified to accomodate registering nodal integer or real fields and
  element-based integer or real fields separately using registerIntNodalField()
  (not yet implemented), registerRealNodalField(), registerIntElementField()
  (not yet implemented), and registerRealElementField(), respectively.
- The Tribol update() function has been modified to take three arguments, the
  last of which is the returned Tribol timestep vote. The new function is:
  update( int cycle, real t, real &dt ). Users should note that in the
  previous API, dt was not modifiable by Tribol. Now it is and serves as the
  Tribol timestep vote. One must be careful not to override their actual dt
  value if this is not desired.
- The simple_tribol.cpp API has a modified Update() function to take a
  modifiable timestep argument. Previously, this function took a non-modifiable
  dt input argument. Be careful to not override the host-code dt with a call to
  this function if that is not desired.
- setPenaltyData() has been deprecated in favor of more explicit routines to set
  the different penalty calculation options. See setKinematicConstantPenalty()
  and setKinematicElementPenalty(). 
- API support for setting different enforcement options has changed. See
  setPenaltyOptions() and setLagrangeMultiplierImplicitOptions().

### Fixed
- Preliminary bugfix for binning of small surface elements.
- Tolerancing in CheckFacePair() based on external testing.
- Unsmooth common-plane + penalty force behavior appears to have been fixed by
  tolerancing.
- Set number of cell nodes and number of cells in registerMesh() prior to
  calling sort routine.
- Changed parameters.gap_tol default in tribol.cpp back to 1.e-10 from 1.e-4 to
  prevent tensile contact solutions in COMMON_PLANE with PENALTY. The
  parameters.gap_tol is no longer used in SINGLE_MORTAR or ALIGNED_MORTAR since
  Lagrange Multiplier contact is determined by the gap AND pressure solution in
  both of those methods. That is, we do not cull the active set based on nodal
  gap alone.

### Code Design
- Design refactors to allow easy to implement features and new methods.
- Use of templated functions for normal and tangential physics application.
- Minimal testing utilities in `src/tribol/utils/TestUtils.hpp`.
- Added simplified API specifically for external/Fortran testing.
- Refactored valid method, model, enforcement, etc., combination checks.
- Refactored the implementation of the Galerkin evaluations and some of the
  finite element functionality.
- Placed all data associated with penalty stiffness calculations on the MeshData
  class.

### Testing
- Added tests to support various aspects of mortar implementations.
- Added COMMON_PLANE + penalty example.
- Added TIED with COMMON_PLANE + penalty test.
- Added some specific tests for computational geometry and mortar weights.
- Added and/or modified tests for element-based penalty and timestep vote
  features.
- Added gap-rate tests and modified tests based on associated refactoring.
