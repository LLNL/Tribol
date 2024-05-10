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
    SLIC_DEBUG("binarySearch: empty array given");
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

  SLIC_DEBUG("binary_search: could not locate value in provided array.");
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
    [=] TRIBOL_HOST_DEVICE (IndexT i) { return array[i] + range[i] < value; },
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

template <typename T, MemorySpace MSPACE>
TRIBOL_HOST_DEVICE void transpose( const ArrayT<T, 2, MSPACE>& in, 
                                   ArrayT<T, 2, MSPACE>& out )
{
  auto h_in = in.shape()[0];
  auto w_in = in.shape()[1];
  auto h_out = out.shape()[0];
  auto w_out = out.shape()[1];

  SLIC_ERROR_IF(h_in != w_out, "Input number of rows does not equal output number of columns.");
  SLIC_ERROR_IF(w_in != h_out, "Input number of columns does not equal output number of rows.");

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