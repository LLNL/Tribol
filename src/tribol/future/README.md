# Tribol Future Contact Library

`tribol::future` is an independent, statically composed contact library. It does
not include or link the legacy Tribol mesh, search, geometry, or physics layers.

## Composition

The public operator combines four policies:

```cpp
using Contact = tribol::future::MfemContactOperator<
    tribol::future::EnergyMortar<tribol::future::NodalConstraints,
                                 tribol::future::QuadraticPenaltyLaw,
                                 tribol::future::EnzymeLinearization>,
    tribol::future::LowOrderRefinedSurface,
    tribol::future::HipExecution,
    tribol::future::BvhSearch>;
```

The concepts in `Concepts.hpp` constrain each layer only at the operations its
caller consumes. Algorithms are duck-typed and do not inherit from a common
base class. `ContactAlgorithmTraits` supplies conservative capability defaults,
so an unrelated custom algorithm does not need mortar-specific declarations.
Common-plane response laws and single-mortar basis policies likewise expose
small operation-based concepts; algorithms do not identify policies by their
concrete C++ types.
Built-in combinations may also be selected through the optional variant factory
in `RuntimeFactory.hpp`; custom combinations use the typed API.

## Lifecycle

Construction fixes the parent mesh and finite-element spaces. The host then
uses the following explicit lifecycle:

```cpp
Contact contact(mesh, coordinates,
                tribol::future::PairedSurfaces{{mortar}, {nonmortar}},
                parameters);
contact.updateInteractions();
contact.updateGeometry(coordinates);
contact.addResidual(residual);
auto derivatives = contact.linearize();
```

- `updateInteractions()` performs search and replaces the frozen interaction
  topology.
- `updateGeometry()` transfers new coordinates while retaining that topology.
- `rebuildGeometry()` rebuilds submesh/LOR/redecomp data and invalidates the
  interactions.
- Residual and derivative evaluation reject missing or stale interactions.
- AMR, repartitioning, and finite-element-space replacement require a new
  contact object.
- Contact operators are intentionally non-copyable and non-movable so returned
  result and linearization views retain stable storage.

## Data Boundaries

`MfemBridge` owns the boundary submesh, optional LOR mesh, redecomp mesh,
transfers, and true-DOF restrictions. `ContactDomain` is the persistent SoA
kernel representation. Its views contain pointers and scalar metadata only.

`ContactRestrictionOperator` implements the primal restriction and its exact
dual transpose. Algorithm `State` vectors are contact-space vectors in the SoA
layout used by `ContactDomain`; host true-DOF values must be transferred through
the bridge restriction before constructing such a state. Element material
fields are indexed by contact-domain element.

## Execution

`SequentialExecution`, `HipExecution`, and `CudaExecution` use the same method
kernels. Device policies use atomic accumulation by default. Wrapping a backend
in `DeterministicExecution<Backend>` uses a precomputed conflict coloring so
independent interactions still execute in parallel while shared-DOF updates are
ordered reproducibly.

Allocation is confined to construction, geometry rebuild, and interaction
rebuild. Named linearization operators and their workspaces are constructed
once and reused. Residual and operator `Mult()` calls reuse persistent
workspaces. Returning total energy requires an explicit device-to-host
synchronization; matrix-free derivative application does not.

## Pressure Staging

For host-defined nodal pressure laws, call `evaluateNodalKinematics()`, compute
pressure on the host or device, and pass either pressure alone for residual-only
evaluation or `ExternalPressureData` for a consistent linearization. No virtual
host callback is invoked by a device kernel.
