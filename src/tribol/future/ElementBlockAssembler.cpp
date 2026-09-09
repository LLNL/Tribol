#include "tribol/future/ElementBlockAssembler.hpp"

#include <stdexcept>

namespace tribol::future {

ElementBlockAssembler::ElementBlockAssembler( const mfem::ParFiniteElementSpace& parent_test_space,
                                              const mfem::ParFiniteElementSpace& parent_trial_space,
                                              const mfem::FiniteElementSpace& contact_test_space,
                                              const mfem::FiniteElementSpace& contact_trial_space )
    : transfer_( parent_test_space, parent_trial_space, contact_test_space, contact_trial_space )
{
}

std::unique_ptr<mfem::HypreParMatrix> ElementBlockAssembler::assemble( const ElementBlockBatch& blocks ) const
{
  if ( blocks.test_elements.size() != blocks.trial_elements.size() ||
       blocks.test_elements.size() != blocks.matrices.size() ) {
    throw std::invalid_argument( "Element-block index and matrix arrays must have equal lengths." );
  }
  auto matrix = transfer_.TransferToParallel( blocks.test_elements, blocks.trial_elements, blocks.matrices, true );
  return std::unique_ptr<mfem::HypreParMatrix>( matrix.release() );
}

}  // namespace tribol::future
