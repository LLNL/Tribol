#ifndef TRIBOL_FUTURE_RESULTS_HPP_
#define TRIBOL_FUTURE_RESULTS_HPP_

#include "tribol/future/Config.hpp"

namespace tribol::future {

struct NodalKinematicsView {
  const mfem::Vector& gap;
  const mfem::Vector& weighted_gap;
  const mfem::Vector& tributary_area;
};

struct QuadratureDiagnosticsView {
  const mfem::Array<int>& pair_offsets;
  const mfem::Vector& points;
  const mfem::Vector& weights;
  const mfem::Vector& gaps;
};

struct ExternalPressureData {
  const mfem::Vector& potential_density;
  const mfem::Vector& pressure;
  const mfem::Vector& pressure_tangent;
};

struct PenaltyLinearization {
  mfem::Operator* dforce_dx{};
  mfem::Operator* dforce_dvelocity{};
};

struct LagrangeMultiplierLinearization {
  mfem::Operator* dforce_dx{};
  mfem::Operator* dforce_dlambda{};
  mfem::Operator* dgap_dx{};
};

struct PenaltyResultView {
  const mfem::Vector& force;
  Real energy{};
};

struct LagrangeMultiplierResultView {
  const mfem::Vector& primal_force;
  const mfem::Vector& constraint_residual;
};

}  // namespace tribol::future

#endif
