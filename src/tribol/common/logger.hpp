// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef COMMON_LOGGER_HPP_
#define COMMON_LOGGER_HPP_

#ifdef USE_SLIC

#include "axom/slic.hpp"  // for logging
#define TRIBOL_ASSERT( x ) SLIC_ASSERT(x)
#define TRIBOL_ERROR(x) SLIC_ERROR(x)
#else

#define TRIBOL_ASSERT( x )
#define TRIBOL_ERROR(x)

#endif

#endif /* COMMON_LOGGER_HPP_ */
