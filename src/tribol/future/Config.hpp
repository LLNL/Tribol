#ifndef TRIBOL_FUTURE_CONFIG_HPP_
#define TRIBOL_FUTURE_CONFIG_HPP_

#include "tribol/config.hpp"

#include "mfem.hpp"

#include <cstdint>

#if defined( __CUDACC__ ) || defined( __HIPCC__ )
#define TRIBOL_FUTURE_HOST_DEVICE __host__ __device__
#else
#define TRIBOL_FUTURE_HOST_DEVICE
#endif

#if defined( __clang__ )
#define TRIBOL_FUTURE_ALWAYS_INLINE __attribute__( ( always_inline ) ) inline
#else
#define TRIBOL_FUTURE_ALWAYS_INLINE inline
#endif

#ifndef TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER
#define TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER 4
#endif

namespace tribol::future {

using Real = mfem::real_t;
using Index = std::int32_t;

inline constexpr int native_high_order_max_order = TRIBOL_FUTURE_NATIVE_HO_MAX_ORDER;

}  // namespace tribol::future

#endif
