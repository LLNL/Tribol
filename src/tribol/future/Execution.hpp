#ifndef TRIBOL_FUTURE_EXECUTION_HPP_
#define TRIBOL_FUTURE_EXECUTION_HPP_

#include "tribol/future/Config.hpp"

#include <utility>

namespace tribol::future {

struct SequentialExecution {
  static constexpr bool uses_device = false;
  static constexpr bool uses_hip = false;
  static constexpr bool uses_cuda = false;
  static constexpr bool deterministic = false;

  template <typename Body>
  static void forAll( int count, Body&& body )
  {
    mfem::forall_switch( false, count, std::forward<Body>( body ) );
  }

  static void synchronize() {}
};

struct HipExecution {
#ifdef MFEM_USE_HIP
  static constexpr bool available = true;
#else
  static constexpr bool available = false;
#endif
  static constexpr bool uses_device = true;
  static constexpr bool uses_hip = true;
  static constexpr bool uses_cuda = false;
  static constexpr bool deterministic = false;

  template <typename Body>
  static void forAll( int count, Body&& body )
  {
    static_assert( available, "HipExecution requires an MFEM HIP build." );
    mfem::forall_switch( true, count, std::forward<Body>( body ) );
  }

  static void synchronize() { MFEM_DEVICE_SYNC; }
};

struct CudaExecution {
#ifdef MFEM_USE_CUDA
  static constexpr bool available = true;
#else
  static constexpr bool available = false;
#endif
  static constexpr bool uses_device = true;
  static constexpr bool uses_hip = false;
  static constexpr bool uses_cuda = true;
  static constexpr bool deterministic = false;

  template <typename Body>
  static void forAll( int count, Body&& body )
  {
    static_assert( available, "CudaExecution requires an MFEM CUDA build." );
    mfem::forall_switch( true, count, std::forward<Body>( body ) );
  }

  static void synchronize() { MFEM_DEVICE_SYNC; }
};

template <typename Backend>
struct DeterministicExecution {
  using backend_type = Backend;

  static constexpr bool uses_device = Backend::uses_device;
  static constexpr bool uses_hip = Backend::uses_hip;
  static constexpr bool uses_cuda = Backend::uses_cuda;
  static constexpr bool deterministic = true;
  static constexpr bool available = [] {
    if constexpr ( requires { Backend::available; } ) {
      return Backend::available;
    } else {
      return true;
    }
  }();

  template <typename Body>
  static void forAll( int count, Body&& body )
  {
    Backend::forAll( count, std::forward<Body>( body ) );
  }

  static void synchronize() { Backend::synchronize(); }
};

using DeterministicSequentialExecution = DeterministicExecution<SequentialExecution>;
using DeterministicHipExecution = DeterministicExecution<HipExecution>;
using DeterministicCudaExecution = DeterministicExecution<CudaExecution>;

struct AtomicAccumulation {
  template <typename T>
  TRIBOL_FUTURE_HOST_DEVICE static void add( T& destination, T value )
  {
    AtomicAdd( destination, value );
  }
};

}  // namespace tribol::future

#endif
