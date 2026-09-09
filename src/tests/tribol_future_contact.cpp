#include "gtest/gtest.h"

#include "tribol/future/TribolFuture.hpp"

#include <array>
#include <cmath>
#include <memory>

namespace tribol::future {
namespace {

using PenaltyAlgorithm = EnergyMortar<NodalConstraints, QuadraticPenaltyLaw, AnalyticLinearization>;
using PenaltyContact =
    MfemContactOperator<PenaltyAlgorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
using MultiplierAlgorithm = EnergyMortar<NodalConstraints, NodalLagrangeMultiplier, AnalyticLinearization>;
using MultiplierContact =
    MfemContactOperator<MultiplierAlgorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
using ExternalAlgorithm = EnergyMortar<NodalConstraints, ExternalPressureLaw, AnalyticLinearization>;
using ExternalContact = MfemContactOperator<ExternalAlgorithm, LowOrderRefinedSurface, SequentialExecution, BvhSearch>;
using QuadraturePenaltyAlgorithm = EnergyMortar<QuadraturePointConstraints, QuadraticPenaltyLaw>;
using QuadraturePenaltyContact = MfemContactOperator<QuadraturePenaltyAlgorithm, LowOrderRefinedSurface,
                                                     SequentialExecution, CartesianProductSearch>;

static_assert( EnforcementLaw<QuadraticPenaltyLaw> );
static_assert( LinearizationPolicy<AnalyticLinearization> );
static_assert( LinearizationPolicy<EnzymeLinearization> );
static_assert( SurfaceDiscretization<LowOrderRefinedSurface> );
static_assert( SurfaceDiscretization<NativeHighOrderSurface<4>> );
static_assert( ExecutionPolicy<SequentialExecution> );
static_assert( ExecutionPolicy<DeterministicSequentialExecution> );
static_assert( SearchPolicy<CartesianProductSearch> );
static_assert( SearchPolicy<BvhSearch> );
static_assert( CommonPlaneNormalLaw<ConstantKinematicPenalty> );
static_assert( CommonPlaneNormalLaw<MaterialKinematicPenalty> );
static_assert( CommonPlaneNormalLaw<TiedNormalResponse> );
static_assert( CommonPlaneNormalLaw<TiedFullResponse> );
static_assert( CommonPlaneRateLaw<NoRatePenalty> );
static_assert( CommonPlaneRateLaw<ConstantRatePenalty> );
static_assert( CommonPlaneRateLaw<PercentageRatePenalty> );
static_assert( CommonPlaneTangentialLaw<FrictionlessTangentialResponse> );
static_assert( CommonPlaneTangentialLaw<ViscousTangentialResponse> );
static_assert( MortarBasisPolicy<PrimalMortarBasis> );
static_assert( MortarBasisPolicy<DualMortarBasis> );
static_assert( ContactAlgorithm<PenaltyAlgorithm> );
static_assert( ContactAlgorithm<MultiplierAlgorithm> );
static_assert( ContactAlgorithm<CommonPlane<>> );
static_assert( ContactAlgorithm<SingleMortar<>> );
static_assert( std::is_trivial_v<InteractionBatch> );
static_assert( std::is_trivially_copyable_v<InteractionBatch> );

struct InvalidPressureLaw {
  Real pressure( Real ) const { return 0.0; }
};

struct ScalarOnlyPressureLaw {
  Real potential( Real gap ) const { return gap; }
  Real pressure( Real gap ) const { return gap; }
  Real tangent( Real ) const { return 1.0; }
};

struct InvalidExecutionPolicy {
  static constexpr bool uses_device = false;
  static constexpr bool uses_hip = false;
  static constexpr bool uses_cuda = false;
  static constexpr bool deterministic = false;
  static void synchronize() {}
};

struct CustomPressureLaw {
  Real stiffness{ 1.0 };

  template <typename T>
  T potential( const T& gap ) const
  {
    return T{ 0.5 * stiffness } * gap * gap;
  }

  template <typename T>
  T pressure( const T& gap ) const
  {
    return T{ stiffness } * gap;
  }

  template <typename T>
  T tangent( const T& ) const
  {
    return T{ stiffness };
  }
};

struct CustomNormalLaw {
  Real stiffness{ 1.0 };

  static constexpr bool ties_normal = false;
  static constexpr bool ties_tangential = false;
  static constexpr bool needs_reference_field = false;
  static constexpr bool needs_material_fields = false;

  TRIBOL_FUTURE_HOST_DEVICE Real penaltyStiffness( int, int, const Real*, const Real* ) const { return stiffness; }
};

struct CustomRateLaw {
  Real scale{};

  static constexpr bool needs_velocity = true;
  static constexpr bool is_conservative = false;

  TRIBOL_FUTURE_HOST_DEVICE Real rateCoefficient( Real normal_stiffness ) const { return scale * normal_stiffness; }
};

struct CustomTangentialLaw {
  static constexpr bool needs_velocity = false;
  static constexpr bool is_conservative = true;

  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE void addTraction( const std::array<T, 3>&, const T&, const std::array<T, 3>&, int,
                                              std::array<T, 3>& ) const
  {
  }
};

struct CustomMortarBasis {
  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE T value( int, bool, bool, const T&, const T&, const T& primal_value ) const
  {
    return primal_value;
  }
};

static_assert( !EnforcementLaw<InvalidPressureLaw> );
static_assert( !EnforcementLaw<ScalarOnlyPressureLaw> );
static_assert( !ExecutionPolicy<InvalidExecutionPolicy> );
static_assert( EnforcementLaw<CustomPressureLaw> );
static_assert( CommonPlaneNormalLaw<CustomNormalLaw> );
static_assert( CommonPlaneRateLaw<CustomRateLaw> );
static_assert( CommonPlaneTangentialLaw<CustomTangentialLaw> );
static_assert( MortarBasisPolicy<CustomMortarBasis> );
static_assert( ContactAlgorithm<EnergyMortar<NodalConstraints, CustomPressureLaw, AnalyticLinearization>> );
static_assert( ContactAlgorithm<CommonPlane<CustomNormalLaw, CustomRateLaw, CustomTangentialLaw>> );
static_assert( ContactAlgorithm<SingleMortar<CustomMortarBasis>> );

struct InvalidContactAlgorithm {
  struct Parameters {};
  struct State {};
  struct Kinematics {};
  struct InteractionStorage {};
  struct Linearization {};
};

static_assert( !ContactAlgorithm<InvalidContactAlgorithm> );

struct IndependentContactAlgorithm {
  struct Parameters {};
  struct State {};
  struct Kinematics {};
  struct InteractionStorage {};
  struct Linearization {};

  explicit IndependentContactAlgorithm( Parameters = {} ) {}

  void buildInteractions( const ContactDomain&, const CandidatePairs&, InteractionStorage& );
  void prepare( const ContactDomain&, const InteractionStorage&, Kinematics& );
  void evaluateKinematics( const ContactDomain&, const InteractionStorage&, const State&, Kinematics& );
  Real addResidual( const ContactDomain&, const InteractionStorage&, const State&, const Kinematics&, mfem::Vector& );
  void applyCoordinateDerivative( const ContactDomain&, const InteractionStorage&, const State&, const Kinematics&,
                                  const mfem::Vector&, mfem::Vector& );
};

static_assert( ContactAlgorithm<IndependentContactAlgorithm> );
static_assert( OperatorContactAlgorithm<IndependentContactAlgorithm> );
static_assert( !ContactAlgorithmTraits<IndependentContactAlgorithm>::capabilities.needs_dual_space );
static_assert( !ContactAlgorithmTraits<IndependentContactAlgorithm>::has_nodal_kinematics );

class FutureContactTest : public testing::Test {
 protected:
  void SetUp() override
  {
    serial_mesh_ = std::make_unique<mfem::Mesh>(
        mfem::Mesh::MakeCartesian2D( 1, 1, mfem::Element::QUADRILATERAL, true, 1.0, 0.1 ) );
    mesh_ = std::make_unique<mfem::ParMesh>( MPI_COMM_WORLD, *serial_mesh_ );
    mesh_->EnsureNodes();
    coordinates_ = dynamic_cast<mfem::ParGridFunction*>( mesh_->GetNodes() );
    ASSERT_NE( coordinates_, nullptr );
  }

  std::unique_ptr<mfem::Mesh> serial_mesh_;
  std::unique_ptr<mfem::ParMesh> mesh_;
  mfem::ParGridFunction* coordinates_{};
};

mfem::Vector upperSurfaceTranslation( const mfem::ParGridFunction& coordinates, Real x_translation, Real y_translation )
{
  auto* space = coordinates.ParFESpace();
  mfem::ParGridFunction local_direction( space );
  local_direction = 0.0;
  for ( int dof = 0; dof < space->GetNDofs(); ++dof ) {
    const int y_dof = space->DofToVDof( dof, 1 );
    if ( coordinates[y_dof] > 0.05 ) {
      local_direction[space->DofToVDof( dof, 0 )] = x_translation;
      local_direction[y_dof] = y_translation;
    }
  }
  mfem::Vector true_direction( space->GetTrueVSize() );
  local_direction.GetTrueDofs( true_direction );
  return true_direction;
}

TEST_F( FutureContactTest, RestrictionAndTransposeAreAdjoint )
{
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } } );
  const auto& restriction = contact.bridge().coordinateRestriction();
  mfem::Vector parent( restriction.Width() );
  mfem::Vector contact_vector( restriction.Height() );
  for ( int i = 0; i < parent.Size(); ++i ) {
    parent[i] = 0.125 * ( i + 1 );
  }
  for ( int i = 0; i < contact_vector.Size(); ++i ) {
    contact_vector[i] = -0.25 + 0.0625 * i;
  }

  mfem::Vector restricted( restriction.Height() );
  mfem::Vector transposed( restriction.Width() );
  restriction.Mult( parent, restricted );
  restriction.MultTranspose( contact_vector, transposed );

  Real left = restricted * contact_vector;
  Real right = parent * transposed;
  MPI_Allreduce( MPI_IN_PLACE, &left, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &right, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( left, right, 1.0e-12 );
}

TEST_F( FutureContactTest, ContactDomainStoresCoordinatesInStructureOfArraysOrder )
{
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } } );
  const auto domain = contact.bridge().domain().view( false );
  ASSERT_EQ( domain.dimension, 2 );
  for ( int element = 0; element < domain.number_of_elements; ++element ) {
    for ( int node = 0; node < domain.numberOfElementDofs( element ); ++node ) {
      const int dof = domain.elementDof( element, node );
      EXPECT_GE( domain.coordinate( 0, dof ), 0.0 );
      EXPECT_LE( domain.coordinate( 0, dof ), 1.0 );
      EXPECT_GE( domain.coordinate( 1, dof ), 0.0 );
      EXPECT_LE( domain.coordinate( 1, dof ), 0.1 );
    }
  }
}

TEST_F( FutureContactTest, ScalarRestrictionAndTransposeAreAdjoint )
{
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } } );
  const auto& restriction = contact.bridge().scalarRestriction();
  mfem::Vector parent( restriction.Width() );
  mfem::Vector contact_vector( restriction.Height() );
  for ( int i = 0; i < parent.Size(); ++i ) {
    parent[i] = 0.2 * ( i + 1 );
  }
  for ( int i = 0; i < contact_vector.Size(); ++i ) {
    contact_vector[i] = 0.3 - 0.05 * i;
  }

  mfem::Vector restricted( restriction.Height() );
  mfem::Vector transposed( restriction.Width() );
  restriction.Mult( parent, restricted );
  restriction.MultTranspose( contact_vector, transposed );

  Real left = restricted * contact_vector;
  Real right = parent * transposed;
  MPI_Allreduce( MPI_IN_PLACE, &left, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &right, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( left, right, 1.0e-12 );
}

TEST_F( FutureContactTest, HighOrderToLorRestrictionAndTransposeAreAdjoint )
{
  mfem::H1_FECollection collection( 2, mesh_->Dimension() );
  mfem::ParFiniteElementSpace space( mesh_.get(), &collection, 2, mfem::Ordering::byVDIM );
  mfem::ParGridFunction high_order_coordinates( &space );
  mfem::VectorFunctionCoefficient curved_coordinates( 2, []( const mfem::Vector& point, mfem::Vector& value ) {
    value.SetSize( 2 );
    value[0] = point[0];
    value[1] = point[1] + 0.03 * point[0] * ( 1.0 - point[0] );
  } );
  high_order_coordinates.ProjectCoefficient( curved_coordinates );

  PenaltyContact contact( *mesh_, high_order_coordinates, PairedSurfaces{ { 3 }, { 1 } } );
  const auto& restriction = contact.bridge().coordinateRestriction();
  mfem::Vector parent( restriction.Width() );
  mfem::Vector contact_vector( restriction.Height() );
  for ( int index = 0; index < parent.Size(); ++index ) {
    parent[index] = 0.125 * ( index + 1 );
  }
  for ( int index = 0; index < contact_vector.Size(); ++index ) {
    contact_vector[index] = 0.4 - 0.025 * index;
  }
  mfem::Vector restricted( restriction.Height() );
  mfem::Vector transposed( restriction.Width() );
  restriction.Mult( parent, restricted );
  restriction.MultTranspose( contact_vector, transposed );

  Real left = restricted * contact_vector;
  Real right = parent * transposed;
  MPI_Allreduce( MPI_IN_PLACE, &left, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &right, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( left, right, 1.0e-11 );

  const auto& scalar_restriction = contact.bridge().scalarRestriction();
  mfem::Vector parent_scalar( scalar_restriction.Width() );
  mfem::Vector contact_scalar( scalar_restriction.Height() );
  for ( int index = 0; index < parent_scalar.Size(); ++index ) {
    parent_scalar[index] = -0.2 + 0.1 * index;
  }
  for ( int index = 0; index < contact_scalar.Size(); ++index ) {
    contact_scalar[index] = 0.3 - 0.04 * index;
  }
  mfem::Vector restricted_scalar( scalar_restriction.Height() );
  mfem::Vector transposed_scalar( scalar_restriction.Width() );
  scalar_restriction.Mult( parent_scalar, restricted_scalar );
  scalar_restriction.MultTranspose( contact_scalar, transposed_scalar );
  left = restricted_scalar * contact_scalar;
  right = parent_scalar * transposed_scalar;
  MPI_Allreduce( MPI_IN_PLACE, &left, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &right, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( left, right, 1.0e-11 );
}

TEST_F( FutureContactTest, NativeHighOrderRejectsOrderAboveConfiguredMaximum )
{
  mfem::H1_FECollection collection( 2, mesh_->Dimension() );
  mfem::ParFiniteElementSpace space( mesh_.get(), &collection, 2, mfem::Ordering::byVDIM );
  mfem::ParGridFunction high_order_coordinates( &space );
  high_order_coordinates.ProjectGridFunction( *coordinates_ );

  using LimitedContact =
      MfemContactOperator<PenaltyAlgorithm, NativeHighOrderSurface<1>, SequentialExecution, CartesianProductSearch>;
  EXPECT_THROW( LimitedContact( *mesh_, high_order_coordinates, PairedSurfaces{ { 3 }, { 1 } } ),
                std::invalid_argument );
}

TEST_F( FutureContactTest, RuntimeFactoryContainsOnlyBuiltInSpecializations )
{
  BuiltInContactParameters parameters;
  parameters.penalty_stiffness = 3.0;
  parameters.search_expansion = 0.2;
  auto contact = makeBuiltInContact( BuiltInContactMethod::EnergyMortarNodalPenalty, *mesh_, *coordinates_,
                                     PairedSurfaces{ { 3 }, { 1 } }, parameters );
  std::visit( []( auto& method ) { method->updateInteractions(); }, contact );
}

TEST_F( FutureContactTest, InteractionRebuildCreatesTopologyAndOrderBatches )
{
  PenaltyContact::Parameters parameters;
  parameters.search.expansion = 0.2;
  MPI_Comm_size( mesh_->GetComm(), &parameters.bridge.redecomp_ranks );
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  int global_pairs = contact.interactions().pairs.Size();
  MPI_Allreduce( MPI_IN_PLACE, &global_pairs, 1, MPI_INT, MPI_SUM, mesh_->GetComm() );
  EXPECT_EQ( global_pairs, 1 );
  if ( contact.interactions().pairs.Size() == 0 ) {
    EXPECT_EQ( contact.interactions().batches.Size(), 0 );
    ASSERT_EQ( contact.interactions().coloring.offsets.Size(), 1 );
    EXPECT_EQ( contact.interactions().coloring.offsets[0], 0 );
  } else {
    ASSERT_EQ( contact.interactions().batches.Size(), 1 );
    const auto batch = contact.interactions().batches[0];
    EXPECT_EQ( batch.begin, 0 );
    EXPECT_EQ( batch.end, contact.interactions().pairs.Size() );
    EXPECT_EQ( batch.mortar_topology, ElementTopology::Segment );
    EXPECT_EQ( batch.nonmortar_topology, ElementTopology::Segment );
    EXPECT_EQ( batch.mortar_order, 1 );
    EXPECT_EQ( batch.nonmortar_order, 1 );
    ASSERT_EQ( contact.interactions().coloring.offsets.Size(), 2 );
    EXPECT_EQ( contact.interactions().coloring.offsets[0], 0 );
    EXPECT_EQ( contact.interactions().coloring.offsets[1], 1 );
    ASSERT_EQ( contact.interactions().coloring.pair_indices.Size(), 1 );
    EXPECT_EQ( contact.interactions().coloring.pair_indices[0], 0 );
  }
}

TEST_F( FutureContactTest, DeterministicExecutionMatchesDefaultSequentialExecution )
{
  using DeterministicContact = MfemContactOperator<PenaltyAlgorithm, LowOrderRefinedSurface,
                                                   DeterministicSequentialExecution, CartesianProductSearch>;
  PenaltyContact::Parameters default_parameters;
  default_parameters.algorithm.enforcement.stiffness = 3.0;
  default_parameters.search.expansion = 0.2;
  DeterministicContact::Parameters deterministic_parameters;
  deterministic_parameters.algorithm.enforcement.stiffness = 3.0;
  deterministic_parameters.search.expansion = 0.2;
  PenaltyContact default_contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, default_parameters );
  DeterministicContact deterministic_contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } },
                                              deterministic_parameters );
  default_contact.updateInteractions();
  deterministic_contact.updateInteractions();

  mfem::Vector default_residual( coordinates_->ParFESpace()->GetTrueVSize() );
  mfem::Vector deterministic_residual( coordinates_->ParFESpace()->GetTrueVSize() );
  default_residual = 0.0;
  deterministic_residual = 0.0;
  default_contact.addResidual( default_residual );
  deterministic_contact.addResidual( deterministic_residual );
  deterministic_residual -= default_residual;
  EXPECT_EQ( deterministic_residual.Norml2(), 0.0 );
}

TEST_F( FutureContactTest, ElementBlocksAssembleThroughRedecomp )
{
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } } );
  const auto& bridge = contact.bridge();
  const auto& contact_space = bridge.contactCoordinateSpace();
  axom::Array<int> test_elements( contact_space.GetNE(), contact_space.GetNE() );
  axom::Array<int> trial_elements( contact_space.GetNE(), contact_space.GetNE() );
  axom::Array<mfem::DenseMatrix> matrices( contact_space.GetNE(), contact_space.GetNE() );
  for ( int element = 0; element < contact_space.GetNE(); ++element ) {
    test_elements[element] = element;
    trial_elements[element] = element;
    mfem::Array<int> dofs;
    contact_space.GetElementVDofs( element, dofs );
    matrices[element].SetSize( dofs.Size() );
    matrices[element] = 0.0;
    for ( int row = 0; row < dofs.Size(); ++row ) {
      matrices[element]( row, row ) = 1.0;
    }
  }
  ElementBlockAssembler assembler( bridge.surfaceCoordinateSpace(), bridge.surfaceCoordinateSpace(), contact_space,
                                   contact_space );
  const auto matrix = assembler.assemble( { test_elements, trial_elements, matrices } );
  ASSERT_NE( matrix, nullptr );
  EXPECT_EQ( matrix->Height(), bridge.surfaceCoordinateSpace().GetTrueVSize() );
  EXPECT_EQ( matrix->Width(), bridge.surfaceCoordinateSpace().GetTrueVSize() );
}

TEST_F( FutureContactTest, NativeHighOrderEnergyMortarUsesAllSegmentNodes )
{
  mfem::H1_FECollection collection( 2, mesh_->Dimension() );
  mfem::ParFiniteElementSpace space( mesh_.get(), &collection, 2, mfem::Ordering::byVDIM );
  mfem::ParGridFunction high_order_coordinates( &space );
  mfem::VectorFunctionCoefficient curved_coordinates( 2, []( const mfem::Vector& point, mfem::Vector& value ) {
    value.SetSize( 2 );
    value[0] = point[0];
    value[1] = point[1] + 0.02 * point[0] * ( 1.0 - point[0] );
  } );
  high_order_coordinates.ProjectCoefficient( curved_coordinates );

  using NativeContact =
      MfemContactOperator<PenaltyAlgorithm, NativeHighOrderSurface<4>, SequentialExecution, CartesianProductSearch>;
  NativeContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 5.0;
  parameters.search.expansion = 0.2;
  NativeContact contact( *mesh_, high_order_coordinates, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  int maximum_nodes{};
  const auto view = contact.bridge().domain().view( false );
  for ( int element = 0; element < view.number_of_elements; ++element ) {
    maximum_nodes = std::max( maximum_nodes, view.numberOfElementDofs( element ) );
  }
  MPI_Allreduce( MPI_IN_PLACE, &maximum_nodes, 1, MPI_INT, MPI_MAX, mesh_->GetComm() );
  EXPECT_EQ( maximum_nodes, 3 );

  const auto nodal = contact.evaluateNodalKinematics();
  EXPECT_EQ( nodal.gap.Size(), contact.bridge().scalarRestriction().Width() );
  mfem::Vector residual( space.GetTrueVSize() );
  residual = 0.0;
  const auto response = contact.addResidual( residual );
  EXPECT_TRUE( std::isfinite( response.energy ) );
}

TEST_F( FutureContactTest, LinearOnlyAlgorithmsRejectNativeHighOrderElements )
{
  mfem::H1_FECollection collection( 2, mesh_->Dimension() );
  mfem::ParFiniteElementSpace space( mesh_.get(), &collection, 2, mfem::Ordering::byVDIM );
  mfem::ParGridFunction high_order_coordinates( &space );
  mfem::VectorFunctionCoefficient identity( 2, []( const mfem::Vector& point, mfem::Vector& value ) {
    value.SetSize( 2 );
    value = point;
  } );
  high_order_coordinates.ProjectCoefficient( identity );

  using Contact =
      MfemContactOperator<CommonPlane<>, NativeHighOrderSurface<4>, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, high_order_coordinates, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  EXPECT_THROW( contact.updateInteractions(), std::invalid_argument );
}

TEST_F( FutureContactTest, SingleMortarRejectsTwoDimensionalSurfaces )
{
  using Contact =
      MfemContactOperator<SingleMortar<>, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  EXPECT_THROW( contact.updateInteractions(), std::invalid_argument );
}

TEST_F( FutureContactTest, ExplicitLifecycleAndNamedLinearizations )
{
  PenaltyContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 10.0;
  parameters.search.expansion = 0.2;
  PenaltyContact penalty( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );

  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  EXPECT_THROW( penalty.addResidual( residual ), std::logic_error );
  penalty.updateInteractions();
  const auto kinematics = penalty.evaluateNodalKinematics();
  EXPECT_EQ( kinematics.gap.Size(), penalty.bridge().domain().numberOfDofs() );
  const auto response = penalty.addResidual( residual );
  EXPECT_TRUE( std::isfinite( response.energy ) );
  auto tangent = penalty.linearize();
  ASSERT_NE( tangent.dforce_dx, nullptr );
  const auto repeated_tangent = penalty.linearize();
  EXPECT_EQ( tangent.dforce_dx, repeated_tangent.dforce_dx );

  MultiplierContact multiplier( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } } );
  multiplier.updateInteractions();
  mfem::Vector dual( multiplier.bridge().scalarRestriction().Width() );
  dual = 1.0;
  residual = 0.0;
  const auto lm_response = multiplier.addResidual( dual, residual );
  EXPECT_EQ( lm_response.constraint_residual.Size(), dual.Size() );
  auto lm_linearization = multiplier.linearize( dual );
  ASSERT_NE( lm_linearization.dforce_dx, nullptr );
  ASSERT_NE( lm_linearization.dforce_dlambda, nullptr );
  ASSERT_NE( lm_linearization.dgap_dx, nullptr );
  const auto repeated_lm_linearization = multiplier.linearize( dual );
  EXPECT_EQ( lm_linearization.dforce_dx, repeated_lm_linearization.dforce_dx );
  EXPECT_EQ( lm_linearization.dforce_dlambda, repeated_lm_linearization.dforce_dlambda );
  EXPECT_EQ( lm_linearization.dgap_dx, repeated_lm_linearization.dgap_dx );
}

TEST_F( FutureContactTest, GeometryUpdatesPreservePairsAndInvalidateOldLinearizations )
{
  PenaltyContact::Parameters parameters;
  parameters.search.expansion = 0.2;
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();
  auto linearization = contact.linearize();
  mfem::Vector direction( coordinates_->ParFESpace()->GetTrueVSize() );
  mfem::Vector result( direction.Size() );
  direction = 0.0;
  contact.updateGeometry( *coordinates_ );
  EXPECT_THROW( linearization.dforce_dx->Mult( direction, result ), std::logic_error );

  contact.rebuildGeometry( *coordinates_ );
  mfem::Vector residual( direction.Size() );
  residual = 0.0;
  EXPECT_THROW( contact.addResidual( residual ), std::logic_error );
}

TEST_F( FutureContactTest, PenaltyForceConservesLinearMomentumAndTangentIsSymmetric )
{
  PenaltyContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 7.0;
  parameters.search.expansion = 0.2;
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();
  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  contact.addResidual( residual );
  mfem::Vector local_residual( coordinates_->Size() );
  coordinates_->ParFESpace()->GetProlongationMatrix()->Mult( residual, local_residual );
  Real component_sum[2]{};
  for ( int component = 0; component < 2; ++component ) {
    for ( int dof = 0; dof < coordinates_->ParFESpace()->GetNDofs(); ++dof ) {
      component_sum[component] += local_residual[coordinates_->ParFESpace()->DofToVDof( dof, component )];
    }
  }
  MPI_Allreduce( MPI_IN_PLACE, component_sum, 2, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( component_sum[0], 0.0, 1.0e-11 );
  EXPECT_NEAR( component_sum[1], 0.0, 1.0e-11 );

  mfem::Vector first( residual.Size() );
  mfem::Vector second( residual.Size() );
  for ( int dof = 0; dof < residual.Size(); ++dof ) {
    first[dof] = 0.01 * ( dof + 1 );
    second[dof] = -0.03 + 0.005 * dof;
  }
  const auto linearization = contact.linearize();
  mfem::Vector first_action( residual.Size() );
  mfem::Vector second_action( residual.Size() );
  linearization.dforce_dx->Mult( first, first_action );
  linearization.dforce_dx->Mult( second, second_action );
  Real first_product = first * second_action;
  Real second_product = second * first_action;
  MPI_Allreduce( MPI_IN_PLACE, &first_product, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &second_product, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( first_product, second_product, 1.0e-9 );
}

TEST_F( FutureContactTest, SelfContactSuppressesDuplicatesAndConfiguredAdjacency )
{
  using Algorithm = CommonPlane<>;
  using Contact = MfemContactOperator<Algorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, SelfContactSurface{ { 1, 2, 3, 4 }, 1 }, parameters );
  contact.updateInteractions();
  const auto view = contact.bridge().domain().view( false );
  for ( const auto pair : contact.interactions().pairs ) {
    EXPECT_NE( pair.mortar_element, pair.nonmortar_element );
    for ( int first = 0; first < view.numberOfElementDofs( pair.mortar_element ); ++first ) {
      for ( int second = 0; second < view.numberOfElementDofs( pair.nonmortar_element ); ++second ) {
        EXPECT_NE( view.elementDof( pair.mortar_element, first ), view.elementDof( pair.nonmortar_element, second ) );
      }
    }
  }
}

TEST_F( FutureContactTest, PenaltyJvpMatchesCenteredDifferenceWithFrozenPairs )
{
  PenaltyContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 7.0;
  parameters.search.expansion = 0.2;
  PenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  mfem::Vector direction = upperSurfaceTranslation( *coordinates_, 0.0, 0.01 );
  auto linearization = contact.linearize();
  mfem::Vector jvp( direction.Size() );
  linearization.dforce_dx->Mult( direction, jvp );

  mfem::Vector local_direction( coordinates_->Size() );
  coordinates_->ParFESpace()->GetProlongationMatrix()->Mult( direction, local_direction );
  constexpr Real step = 1.0e-6;
  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector plus( direction.Size() );
  plus = 0.0;
  contact.addResidual( plus );

  coordinates_->Add( -2.0 * step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector minus( direction.Size() );
  minus = 0.0;
  contact.addResidual( minus );

  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  plus -= minus;
  plus /= 2.0 * step;
  EXPECT_LT( std::abs( plus.Norml2() - jvp.Norml2() ), 1.0e-6 );
  plus -= jvp;
  EXPECT_LT( plus.Norml2(), 1.0e-6 );
}

TEST_F( FutureContactTest, QuadraturePenaltyDefaultLinearizationProducesFiniteJvp )
{
  QuadraturePenaltyContact::Parameters parameters;
  parameters.algorithm.enforcement.stiffness = 7.0;
  parameters.algorithm.quadrature_points = 3;
  parameters.search.expansion = 0.2;
  QuadraturePenaltyContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();
  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  contact.addResidual( residual );
  auto linearization = contact.linearize();
  mfem::Vector direction = upperSurfaceTranslation( *coordinates_, 0.0, 0.01 );
  mfem::Vector derivative( residual.Size() );
  linearization.dforce_dx->Mult( direction, derivative );
  for ( int dof = 0; dof < derivative.Size(); ++dof ) {
    EXPECT_TRUE( std::isfinite( derivative[dof] ) );
  }

  mfem::Vector local_direction( coordinates_->Size() );
  coordinates_->ParFESpace()->GetProlongationMatrix()->Mult( direction, local_direction );
  constexpr Real step = 1.0e-6;
  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector plus( residual.Size() );
  plus = 0.0;
  contact.addResidual( plus );
  coordinates_->Add( -2.0 * step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector minus( residual.Size() );
  minus = 0.0;
  contact.addResidual( minus );
  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  plus -= minus;
  plus /= 2.0 * step;
  plus -= derivative;
  EXPECT_LT( plus.Norml2(), 1.0e-6 );
}

TEST_F( FutureContactTest, HostPressureUsesExplicitStagedData )
{
  ExternalContact::Parameters parameters;
  parameters.search.expansion = 0.2;
  ExternalContact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();
  const auto nodal = contact.evaluateNodalKinematics();
  mfem::Vector potential( nodal.gap.Size() );
  mfem::Vector pressure( nodal.gap.Size() );
  mfem::Vector tangent( nodal.gap.Size() );
  for ( int dof = 0; dof < nodal.gap.Size(); ++dof ) {
    const Real active_gap = std::min( nodal.gap[dof], 0.0 );
    potential[dof] = 5.0 * active_gap * active_gap;
    pressure[dof] = 10.0 * active_gap;
    tangent[dof] = nodal.gap[dof] < 0.0 ? 10.0 : 0.0;
  }
  const ExternalPressureData data{ potential, pressure, tangent };
  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  const auto result = contact.addResidual( pressure, residual );
  EXPECT_TRUE( std::isfinite( result.energy ) );
  const auto linearization = contact.linearize( data );
  ASSERT_NE( linearization.dforce_dx, nullptr );
}

TEST_F( FutureContactTest, CommonPlaneUsesIndependentAlgorithmInterface )
{
  using Algorithm =
      CommonPlane<ConstantKinematicPenalty, NoRatePenalty, FrictionlessTangentialResponse, AnalyticLinearization>;
  using Contact = MfemContactOperator<Algorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.algorithm.normal.stiffness = 2.0;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();
  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  Algorithm::State state;
  const auto result = contact.addResidual( state, residual );
  EXPECT_TRUE( std::isfinite( result.energy ) );
  auto linearization = contact.linearize( state );
  ASSERT_NE( linearization.dforce_dx, nullptr );
}

TEST_F( FutureContactTest, CommonPlaneUsesProjectedOverlapAfterLargeTangentialMotion )
{
  using Algorithm = CommonPlane<>;
  using Contact = MfemContactOperator<Algorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  const mfem::Vector translation = upperSurfaceTranslation( *coordinates_, 0.25, 0.0 );
  mfem::Vector local_translation( coordinates_->Size() );
  coordinates_->ParFESpace()->GetProlongationMatrix()->Mult( translation, local_translation );
  *coordinates_ += local_translation;
  contact.updateGeometry( *coordinates_ );

  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  contact.addResidual( Algorithm::State{}, residual );
  const auto diagnostics = contact.quadratureDiagnostics();
  Real overlap{};
  for ( int point = 0; point < diagnostics.weights.Size(); ++point ) {
    overlap += diagnostics.weights[point];
  }
  MPI_Allreduce( MPI_IN_PLACE, &overlap, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_NEAR( overlap, 0.75, 1.0e-12 );
}

TEST_F( FutureContactTest, CommonPlaneTiedPoliciesUseReferenceDisplacement )
{
  using TiedNormalAlgorithm =
      CommonPlane<TiedNormalResponse, NoRatePenalty, FrictionlessTangentialResponse, AnalyticLinearization>;
  using TiedFullAlgorithm =
      CommonPlane<TiedFullResponse, NoRatePenalty, FrictionlessTangentialResponse, AnalyticLinearization>;
  using TiedNormalContact =
      MfemContactOperator<TiedNormalAlgorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  using TiedFullContact =
      MfemContactOperator<TiedFullAlgorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;

  TiedNormalContact::Parameters normal_parameters;
  normal_parameters.search.expansion = 0.2;
  TiedFullContact::Parameters full_parameters;
  full_parameters.search.expansion = 0.2;
  TiedNormalContact tied_normal( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, normal_parameters );
  TiedFullContact tied_full( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, full_parameters );
  tied_normal.updateInteractions();
  tied_full.updateInteractions();

  mfem::Vector reference( tied_normal.bridge().domain().coordinateSize() );
  reference = tied_normal.bridge().domain().coordinates();
  const auto domain = tied_normal.bridge().domain().view( false );
  int global_pairs = tied_normal.interactions().pairs.Size();
  MPI_Allreduce( MPI_IN_PLACE, &global_pairs, 1, MPI_INT, MPI_SUM, mesh_->GetComm() );
  ASSERT_EQ( global_pairs, 1 );
  if ( tied_normal.interactions().pairs.Size() > 0 ) {
    const auto pair = tied_normal.interactions().pairs[0];
    const auto surface_normal = [&]( int element ) {
      const int first_dof = domain.elementDof( element, 0 );
      const int second_dof = domain.elementDof( element, 1 );
      const Real dx = domain.coordinate( 0, second_dof ) - domain.coordinate( 0, first_dof );
      const Real dy = domain.coordinate( 1, second_dof ) - domain.coordinate( 1, first_dof );
      const Real length = std::sqrt( dx * dx + dy * dy );
      return std::array<Real, 2>{ dy / length, -dx / length };
    };
    const auto mortar_normal = surface_normal( pair.mortar_element );
    const auto nonmortar_normal = surface_normal( pair.nonmortar_element );
    std::array<Real, 2> common_normal{ nonmortar_normal[0] - mortar_normal[0], nonmortar_normal[1] - mortar_normal[1] };
    Real common_normal_norm = std::sqrt( common_normal[0] * common_normal[0] + common_normal[1] * common_normal[1] );
    if ( common_normal_norm < 1.0e-14 ) {
      common_normal = nonmortar_normal;
      common_normal_norm = 1.0;
    }
    const std::array<Real, 2> tangent{ -common_normal[1] / common_normal_norm, common_normal[0] / common_normal_norm };
    for ( int element = 0; element < domain.number_of_elements; ++element ) {
      if ( domain.side( element ) != SurfaceSide::Mortar ) {
        continue;
      }
      for ( int node = 0; node < domain.numberOfElementDofs( element ); ++node ) {
        const int dof = domain.elementDof( element, node );
        reference[dof] -= 0.05 * tangent[0];
        reference[domain.number_of_dofs + dof] -= 0.05 * tangent[1];
      }
    }
  }

  TiedNormalAlgorithm::State normal_state;
  mfem::Vector normal_residual( coordinates_->ParFESpace()->GetTrueVSize() );
  normal_residual = 0.0;
  EXPECT_THROW( tied_normal.addResidual( normal_state, normal_residual ), std::invalid_argument );
  normal_state.reference_coordinates = &reference;
  tied_normal.addResidual( normal_state, normal_residual );

  TiedFullAlgorithm::State full_state;
  full_state.reference_coordinates = &reference;
  mfem::Vector full_residual( coordinates_->ParFESpace()->GetTrueVSize() );
  full_residual = 0.0;
  tied_full.addResidual( full_state, full_residual );

  Real normal_norm_squared = normal_residual * normal_residual;
  Real full_norm_squared = full_residual * full_residual;
  MPI_Allreduce( MPI_IN_PLACE, &normal_norm_squared, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  MPI_Allreduce( MPI_IN_PLACE, &full_norm_squared, 1, MPI_DOUBLE, MPI_SUM, mesh_->GetComm() );
  EXPECT_LT( std::sqrt( normal_norm_squared ), 1.0e-12 );
  EXPECT_GT( std::sqrt( full_norm_squared ), 1.0e-8 );
}

TEST_F( FutureContactTest, CommonPlaneMaterialFieldsAreValidatedIndependently )
{
  using Algorithm =
      CommonPlane<MaterialKinematicPenalty, NoRatePenalty, FrictionlessTangentialResponse, AnalyticLinearization>;
  using Contact = MfemContactOperator<Algorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  Algorithm::State state;
  mfem::Vector residual( coordinates_->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  EXPECT_THROW( contact.addResidual( state, residual ), std::invalid_argument );

  mfem::Vector thickness( contact.bridge().domain().numberOfElements() );
  mfem::Vector modulus( contact.bridge().domain().numberOfElements() );
  thickness = 0.1;
  modulus = 2.0;
  state.element_thickness = &thickness;
  state.material_modulus = &modulus;
  EXPECT_NO_THROW( contact.addResidual( state, residual ) );
}

TEST_F( FutureContactTest, CommonPlaneRateAndViscousDerivativesMatchFiniteDifferences )
{
  using Algorithm =
      CommonPlane<ConstantKinematicPenalty, ConstantRatePenalty, ViscousTangentialResponse, AnalyticLinearization>;
  using Contact = MfemContactOperator<Algorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.algorithm.normal.stiffness = 4.0;
  parameters.algorithm.rate.coefficient = 0.7;
  parameters.algorithm.tangential.damping = 0.3;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh_, *coordinates_, PairedSurfaces{ { 3 }, { 1 } }, parameters );
  contact.updateInteractions();

  mfem::Vector coordinate_direction = upperSurfaceTranslation( *coordinates_, 0.0, 0.01 );
  mfem::Vector velocity_direction = upperSurfaceTranslation( *coordinates_, 0.01, 0.0 );
  mfem::Vector contact_velocity_direction( contact.bridge().domain().coordinateSize() );
  contact.bridge().coordinateRestriction().Mult( velocity_direction, contact_velocity_direction );
  mfem::Vector velocity( contact.bridge().domain().coordinateSize() );
  contact.bridge().coordinateRestriction().Mult( upperSurfaceTranslation( *coordinates_, 0.02, 0.03 ), velocity );
  Algorithm::State state;
  state.velocity = &velocity;

  const auto linearization = contact.linearize( state );
  mfem::Vector coordinate_jvp( coordinate_direction.Size() );
  mfem::Vector velocity_jvp( velocity_direction.Size() );
  linearization.dforce_dx->Mult( coordinate_direction, coordinate_jvp );
  linearization.dforce_dvelocity->Mult( velocity_direction, velocity_jvp );

  constexpr Real step = 1.0e-6;
  mfem::Vector local_direction( coordinates_->Size() );
  coordinates_->ParFESpace()->GetProlongationMatrix()->Mult( coordinate_direction, local_direction );
  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector coordinate_plus( coordinate_direction.Size() );
  coordinate_plus = 0.0;
  contact.addResidual( state, coordinate_plus );
  coordinates_->Add( -2.0 * step, local_direction );
  contact.updateGeometry( *coordinates_ );
  mfem::Vector coordinate_minus( coordinate_direction.Size() );
  coordinate_minus = 0.0;
  contact.addResidual( state, coordinate_minus );
  coordinates_->Add( step, local_direction );
  contact.updateGeometry( *coordinates_ );
  coordinate_plus -= coordinate_minus;
  coordinate_plus /= 2.0 * step;
  coordinate_plus -= coordinate_jvp;
  EXPECT_LT( coordinate_plus.Norml2(), 1.0e-6 );

  velocity.Add( step, contact_velocity_direction );
  mfem::Vector velocity_plus( velocity_direction.Size() );
  velocity_plus = 0.0;
  contact.addResidual( state, velocity_plus );
  velocity.Add( -2.0 * step, contact_velocity_direction );
  mfem::Vector velocity_minus( velocity_direction.Size() );
  velocity_minus = 0.0;
  contact.addResidual( state, velocity_minus );
  velocity.Add( step, contact_velocity_direction );
  velocity_plus -= velocity_minus;
  velocity_plus /= 2.0 * step;
  velocity_plus -= velocity_jvp;
  EXPECT_LT( velocity_plus.Norml2(), 1.0e-6 );
}

TEST( FutureContactThreeDimensionalTest, CommonPlaneAndSingleMortarAcceptLinearFaces )
{
  auto serial_mesh =
      std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian3D( 1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 0.1 ) );
  auto mesh = std::make_unique<mfem::ParMesh>( MPI_COMM_WORLD, *serial_mesh );
  mesh->EnsureNodes();
  auto* coordinates = dynamic_cast<mfem::ParGridFunction*>( mesh->GetNodes() );
  ASSERT_NE( coordinates, nullptr );

  using CpAlgorithm =
      CommonPlane<ConstantKinematicPenalty, ConstantRatePenalty, ViscousTangentialResponse, AnalyticLinearization>;
  using CpContact =
      MfemContactOperator<CpAlgorithm, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  CpContact::Parameters cp_parameters;
  cp_parameters.algorithm.normal.stiffness = 4.0;
  cp_parameters.algorithm.rate.coefficient = 0.2;
  cp_parameters.algorithm.tangential.damping = 0.1;
  cp_parameters.search.expansion = 0.2;
  CpContact common_plane( *mesh, *coordinates, PairedSurfaces{ { 6 }, { 5 } }, cp_parameters );
  common_plane.updateInteractions();
  CpAlgorithm::State cp_state;
  mfem::Vector contact_velocity( common_plane.bridge().domain().coordinateSize() );
  contact_velocity = 0.0;
  cp_state.velocity = &contact_velocity;
  mfem::Vector residual( coordinates->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  common_plane.addResidual( cp_state, residual );
  auto cp_linearization = common_plane.linearize( cp_state );
  ASSERT_NE( cp_linearization.dforce_dx, nullptr );
  ASSERT_NE( cp_linearization.dforce_dvelocity, nullptr );

  using Mortar =
      MfemContactOperator<SingleMortar<>, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Mortar::Parameters mortar_parameters;
  mortar_parameters.search.expansion = 0.2;
  Mortar mortar( *mesh, *coordinates, PairedSurfaces{ { 6 }, { 5 } }, mortar_parameters );
  mortar.updateInteractions();
  mfem::Vector multiplier( mortar.bridge().scalarRestriction().Width() );
  multiplier = 1.0;
  residual = 0.0;
  const auto result = mortar.addResidual( multiplier, residual );
  EXPECT_EQ( result.constraint_residual.Size(), multiplier.Size() );
  const auto blocks = mortar.linearize( multiplier );
  ASSERT_NE( blocks.dforce_dx, nullptr );
  ASSERT_NE( blocks.dforce_dlambda, nullptr );
  ASSERT_NE( blocks.dgap_dx, nullptr );
}

TEST( FutureContactThreeDimensionalTest, CommonPlaneAndSingleMortarAcceptTriangularFaces )
{
  auto serial_mesh =
      std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian3D( 1, 1, 1, mfem::Element::TETRAHEDRON, 1.0, 1.0, 0.1 ) );
  auto mesh = std::make_unique<mfem::ParMesh>( MPI_COMM_WORLD, *serial_mesh );
  mesh->EnsureNodes();
  auto* coordinates = dynamic_cast<mfem::ParGridFunction*>( mesh->GetNodes() );
  ASSERT_NE( coordinates, nullptr );

  using CommonPlaneContact =
      MfemContactOperator<CommonPlane<>, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  CommonPlaneContact::Parameters common_plane_parameters;
  common_plane_parameters.search.expansion = 0.2;
  CommonPlaneContact common_plane( *mesh, *coordinates, PairedSurfaces{ { 6 }, { 1 } }, common_plane_parameters );
  common_plane.updateInteractions();
  int common_plane_pairs = common_plane.interactions().pairs.Size();
  MPI_Allreduce( MPI_IN_PLACE, &common_plane_pairs, 1, MPI_INT, MPI_SUM, mesh->GetComm() );
  EXPECT_GT( common_plane_pairs, 0 );
  mfem::Vector residual( coordinates->ParFESpace()->GetTrueVSize() );
  residual = 0.0;
  EXPECT_NO_THROW( common_plane.addResidual( CommonPlane<>::State{}, residual ) );

  using MortarContact =
      MfemContactOperator<SingleMortar<>, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  MortarContact::Parameters mortar_parameters;
  mortar_parameters.search.expansion = 0.2;
  MortarContact mortar( *mesh, *coordinates, PairedSurfaces{ { 6 }, { 1 } }, mortar_parameters );
  mortar.updateInteractions();
  int mortar_pairs = mortar.interactions().pairs.Size();
  MPI_Allreduce( MPI_IN_PLACE, &mortar_pairs, 1, MPI_INT, MPI_SUM, mesh->GetComm() );
  EXPECT_GT( mortar_pairs, 0 );
  mfem::Vector dual( mortar.bridge().scalarRestriction().Width() );
  dual = 1.0;
  residual = 0.0;
  EXPECT_NO_THROW( mortar.addResidual( dual, residual ) );
}

TEST( FutureContactThreeDimensionalTest, SingleMortarDerivativeBlocksMatchFiniteDifferences )
{
  auto serial_mesh =
      std::make_unique<mfem::Mesh>( mfem::Mesh::MakeCartesian3D( 1, 1, 1, mfem::Element::HEXAHEDRON, 1.0, 1.0, 0.1 ) );
  auto mesh = std::make_unique<mfem::ParMesh>( MPI_COMM_WORLD, *serial_mesh );
  mesh->EnsureNodes();
  auto* coordinates = dynamic_cast<mfem::ParGridFunction*>( mesh->GetNodes() );
  ASSERT_NE( coordinates, nullptr );

  using Contact =
      MfemContactOperator<SingleMortar<>, LowOrderRefinedSurface, SequentialExecution, CartesianProductSearch>;
  Contact::Parameters parameters;
  parameters.search.expansion = 0.2;
  Contact contact( *mesh, *coordinates, PairedSurfaces{ { 6 }, { 5 } }, parameters );
  contact.updateInteractions();

  mfem::Vector multiplier( contact.bridge().scalarRestriction().Width() );
  mfem::Vector multiplier_direction( multiplier.Size() );
  for ( int dof = 0; dof < multiplier.Size(); ++dof ) {
    multiplier[dof] = 1.0 + 0.1 * dof;
    multiplier_direction[dof] = -0.2 + 0.03 * dof;
  }
  mfem::Vector coordinate_direction( coordinates->ParFESpace()->GetTrueVSize() );
  for ( int dof = 0; dof < coordinate_direction.Size(); ++dof ) {
    coordinate_direction[dof] = 0.005 * ( dof + 1 );
  }

  const auto blocks = contact.linearize( multiplier );
  mfem::Vector force_coordinate_jvp( coordinate_direction.Size() );
  mfem::Vector force_multiplier_jvp( coordinate_direction.Size() );
  mfem::Vector gap_coordinate_jvp( multiplier.Size() );
  blocks.dforce_dx->Mult( coordinate_direction, force_coordinate_jvp );
  blocks.dforce_dlambda->Mult( multiplier_direction, force_multiplier_jvp );
  blocks.dgap_dx->Mult( coordinate_direction, gap_coordinate_jvp );

  constexpr Real step = 1.0e-6;
  mfem::Vector local_direction( coordinates->Size() );
  coordinates->ParFESpace()->GetProlongationMatrix()->Mult( coordinate_direction, local_direction );
  coordinates->Add( step, local_direction );
  contact.updateGeometry( *coordinates );
  mfem::Vector force_plus( coordinate_direction.Size() );
  force_plus = 0.0;
  const auto plus_result = contact.addResidual( multiplier, force_plus );
  mfem::Vector gap_plus( plus_result.constraint_residual );
  coordinates->Add( -2.0 * step, local_direction );
  contact.updateGeometry( *coordinates );
  mfem::Vector force_minus( coordinate_direction.Size() );
  force_minus = 0.0;
  const auto minus_result = contact.addResidual( multiplier, force_minus );
  mfem::Vector gap_minus( minus_result.constraint_residual );
  coordinates->Add( step, local_direction );
  contact.updateGeometry( *coordinates );
  force_plus -= force_minus;
  force_plus /= 2.0 * step;
  force_plus -= force_coordinate_jvp;
  EXPECT_LT( force_plus.Norml2(), 1.0e-6 );
  gap_plus -= gap_minus;
  gap_plus /= 2.0 * step;
  gap_plus -= gap_coordinate_jvp;
  EXPECT_LT( gap_plus.Norml2(), 1.0e-6 );

  mfem::Vector multiplier_plus( multiplier );
  mfem::Vector multiplier_minus( multiplier );
  multiplier_plus.Add( step, multiplier_direction );
  multiplier_minus.Add( -step, multiplier_direction );
  force_plus = 0.0;
  force_minus = 0.0;
  contact.addResidual( multiplier_plus, force_plus );
  contact.addResidual( multiplier_minus, force_minus );
  force_plus -= force_minus;
  force_plus /= 2.0 * step;
  force_plus -= force_multiplier_jvp;
  EXPECT_LT( force_plus.Norml2(), 1.0e-8 );
}

}  // namespace
}  // namespace tribol::future

int main( int argc, char* argv[] )
{
  MPI_Init( &argc, &argv );
  ::testing::InitGoogleTest( &argc, argv );
  const int result = RUN_ALL_TESTS();
  MPI_Finalize();
  return result;
}
