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

  /*!
   * \brief Get instance 
   *
   * \return DataManager instance
   */
  static DataManager& getInstance()
  {
    if (instance_ == nullptr)
    {
      instance_ = std::make_unique<DataManager>();
    }
    return *instance_;
  }

  /*!
   * \brief Get data
   *
   * \param [in] id Data identifier
   *
   * \return Data
   */
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

  /*!
   * \brief Deletes existing data with given id and adds new data
   *
   * \param [in] meshID New mesh id
   *
   * \return MeshData object
   */
  T& addData(IndexT id, T&& data)
  {
    data_map_.erase(id);
    auto data_it = data_map_.emplace(id, std::move(data));
    return *data_it.first;
  }

  /*!
   * \brief Returns pointer to data if found, nullptr otherwise
   *
   * \param [in] id Id of mesh
   *
   * \return True if mesh exists
   */
  T* findData(IndexT id) const
  {
    auto data_it = data_map_.find(id);
    if (data_it == data_map_.end())
    {
      return nullptr;
    } 
    else
    {
      return &(*data_it);
    }
  }

private:
  DataManager() = default;
  ~DataManager() = default;

  std::unordered_map<IndexT, T> data_map_; ///< Unordered map of mesh instances

  static std::unique_ptr<DataManager> instance_;

};

} // end namespace tribol

#endif /* SRC_UTILS_DATAMANAGER_HPP_ */
