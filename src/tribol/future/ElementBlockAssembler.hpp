#ifndef TRIBOL_FUTURE_ELEMENTBLOCKASSEMBLER_HPP_
#define TRIBOL_FUTURE_ELEMENTBLOCKASSEMBLER_HPP_

#include "tribol/future/Config.hpp"

#include "axom/core/Array.hpp"
#include "mfem.hpp"
#include "redecomp/transfer/MatrixTransfer.hpp"

#include <memory>

namespace tribol::future {

struct ElementBlockBatch {
  const axom::Array<int>& test_elements;
  const axom::Array<int>& trial_elements;
  const axom::Array<mfem::DenseMatrix>& matrices;
};

class ElementBlockAssembler {
 public:
  ElementBlockAssembler( const mfem::ParFiniteElementSpace& parent_test_space,
                         const mfem::ParFiniteElementSpace& parent_trial_space,
                         const mfem::FiniteElementSpace& contact_test_space,
                         const mfem::FiniteElementSpace& contact_trial_space );

  std::unique_ptr<mfem::HypreParMatrix> assemble( const ElementBlockBatch& blocks ) const;

 private:
  redecomp::MatrixTransfer transfer_;
};

}  // namespace tribol::future

#endif
