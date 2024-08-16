// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_UTILS_ALGORITHM_HPP_
#define SRC_UTILS_ALGORITHM_HPP_

#include "tribol/common/ArrayTypes.hpp"

namespace tribol
{

namespace algorithm
{

template <typename LCOMP, typename HCOMP>
TRIBOL_HOST_DEVICE IndexT binarySearch( IndexT size,
                                        LCOMP&& lo_comparison,
                                        HCOMP&& hi_comparison )
{
  if (size == 0)
  {
#ifdef TRIBOL_USE_HOST
    SLIC_DEBUG("binarySearch: empty array given");
#endif
    return -1;
  }

  IndexT l = 0;
  IndexT r = size - 1;
  while (l <= r)
  {
    IndexT m = (l + r) / 2;
    if (lo_comparison(m))
    {
      l = m + 1;
    }
    else if (hi_comparison(m))
    {
      r = m - 1;
    }
    else
    {
      return m;
    }
  }

#ifdef TRIBOL_USE_HOST
  SLIC_DEBUG("binary_search: could not locate value in provided array.");
#endif
  return -1;
}

template <typename T>
TRIBOL_HOST_DEVICE IndexT binarySearch( const T* array,
                                        const T* range,
                                        IndexT size,
                                        IndexT value )
{
  return binarySearch(
    size, 
    [=] TRIBOL_HOST_DEVICE (IndexT i) { return array[i] + range[i] <= value; },
    [=] TRIBOL_HOST_DEVICE (IndexT i) { return array[i] > value; }
  );
}

template <typename ARRAY>
TRIBOL_HOST_DEVICE IndexT binarySearch( const ARRAY& array,
                                        const ARRAY& range,
                                        IndexT value )
{
  return binarySearch(array.data(), range.data(), array.size(), value);
}

TRIBOL_HOST_DEVICE inline IndexT symmMatrixRow( IndexT value,
                                                IndexT matrix_width )
{
  return binarySearch(
    matrix_width,
    [=] TRIBOL_HOST_DEVICE (IndexT i) { return (i + 1) * (i + 2) / 2 <= value; },
    [=] TRIBOL_HOST_DEVICE (IndexT i) { return i * (i + 1) / 2 > value; });
}

template <MemorySpace MSPACE, typename T>
TRIBOL_HOST_DEVICE void transpose( const ArrayT<T, 2, MSPACE>& in, 
                                   ArrayT<T, 2, MSPACE>& out )
{
  auto h_in = in.shape()[0];
  auto w_in = in.shape()[1];

#ifdef TRIBOL_USE_HOST
  SLIC_ERROR_IF(h_in != out.shape()[1], "Input number of rows does not equal output number of columns.");
  SLIC_ERROR_IF(w_in != out.shape()[0], "Input number of columns does not equal output number of rows.");
#endif

  for (IndexT i{0}; i < h_in; ++i)
  {
    for (IndexT j{0}; j < w_in; ++j)
    {
      out(j, i) = in(i, j);
    }
  }
}

} // namespace algorithm

} // namespace tribol


#endif /* SRC_UTILS_ALGORITHM_HPP_ */