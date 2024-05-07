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

} // namespace algorithm

} // namespace tribol


#endif /* SRC_UTILS_ALGORITHM_HPP_ */