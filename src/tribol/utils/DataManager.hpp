// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_UTILS_DATAMANAGER_HPP_
#define SRC_UTILS_DATAMANAGER_HPP_

// C/C++ includes
#include <memory>
#include <unordered_map>

// Tribol includes
#include "tribol/types.hpp"

// Axom includes
#include "axom/fmt.hpp"
#include "axom/slic.hpp"

namespace tribol
{

template <typename T>
class DataManager
{
public:
  DataManager(const DataManager& other) = delete;
  DataManager(DataManager&& other) = delete;
  void operator=(const DataManager& other) = delete;
  void operator=(DataManager&& other) = delete;

  static DataManager& getInstance()
  {
    static DataManager instance_;
    return instance_;
  }

  T& at(IndexT id)
  {
    return data_map_.at(id);
  }

  typename std::unordered_map<IndexT, T>::iterator begin() noexcept
  {
    return data_map_.begin();
  }

  typename std::unordered_map<IndexT, T>::iterator end() noexcept
  {
    return data_map_.end();
  }

  size_t erase(IndexT id)
  {
    return data_map_.erase(id);
  }

  void clear() noexcept
  {
    data_map_.clear();
  }

  size_t size() noexcept
  {
    return data_map_.size();
  }

  T& addData(IndexT id, T&& data)
  {
    data_map_.erase(id);
    auto data_it = data_map_.emplace(id, std::move(data));
    return data_it.first->second;
  }

  T* findData(IndexT id)
  {
    auto data_it = data_map_.find(id);
    if (data_it == data_map_.end())
    {
      return nullptr;
    } 
    else
    {
      return &data_it->second;
    }
  }

  T& getData(IndexT id)
  {
    auto data_it = data_map_.find(id);
    SLIC_ERROR_ROOT_IF(data_it == data_map_.end(),
      axom::fmt::format("No data exists for id = {}.", id));
    return data_it->second;
  }

private:
  DataManager() = default;
  ~DataManager() = default;

  std::unordered_map<IndexT, T> data_map_;
};

} // end namespace tribol

#endif /* SRC_UTILS_DATAMANAGER_HPP_ */
