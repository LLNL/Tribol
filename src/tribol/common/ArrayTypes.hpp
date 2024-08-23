// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef TRIBOL_COMMON_ARRAYTYPES_HPP_
#define TRIBOL_COMMON_ARRAYTYPES_HPP_

// Tribol includes
#include "tribol/common/ExecModel.hpp"

#include "axom/core/Array.hpp"
#include "axom/core/ArrayView.hpp"

namespace tribol
{

template <typename T, int DIM, typename ArrayType>
using ArrayBaseT = axom::ArrayBase<T, DIM, ArrayType>;

template <typename T, int DIM = 1, MemorySpace SPACE = MemorySpace::Dynamic>
using ArrayT = axom::Array<T, DIM, toAxomMemorySpace<SPACE>::value>;

template <typename T, int DIM = 1, MemorySpace SPACE = MemorySpace::Dynamic>
using ArrayViewT = axom::ArrayView<T, DIM, toAxomMemorySpace<SPACE>::value>;

template <typename T, int N>
using StackArrayT = axom::StackArray<T, N>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using Array1D = ArrayT<T, 1, SPACE>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using Array2D = ArrayT<T, 2, SPACE>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using Array1DView = ArrayViewT<T, 1, SPACE>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using Array2DView = ArrayViewT<T, 2, SPACE>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using MultiArrayView = ArrayT<ArrayViewT<T, 1, SPACE>, 1, SPACE>;

template <typename T, MemorySpace SPACE = MemorySpace::Dynamic>
using MultiViewArrayView = ArrayViewT<const ArrayViewT<T, 1, SPACE>, 1, SPACE>;

} // namespace tribol

#endif /* TRIBOL_COMMON_ARRAYTYPES_HPP_ */
