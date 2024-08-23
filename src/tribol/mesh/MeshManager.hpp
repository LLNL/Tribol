// Copyright (c) 2017-2023, Lawrence Livermore National Security, LLC and
// other Tribol Project Developers. See the top-level LICENSE file for details.
//
// SPDX-License-Identifier: (MIT)

#ifndef SRC_MESH_MESHMANAGER_HPP_
#define SRC_MESH_MESHMANAGER_HPP_

// Tribol includes
#include "tribol/utils/DataManager.hpp"
#include "tribol/mesh/MeshData.hpp"

namespace tribol
{

// template <MemorySpace MSPACE>


// class MeshManager
// {
// public:
//   MeshManager(const MeshManager& other) = delete;
//   MeshManager(MeshManager&& other) = delete;
//   void operator=(const MeshManager& other) = delete;
//   void operator=(MeshManager&& other) = delete;

//   static MeshManager& getInstance()
//   {
//     static MeshManager instance_;
//     return instance_;
//   }

//   template <MemorySpace MSPACE>
//   MeshData<MSPACE>& at(IndexT id)
//   {
//     return DataManager<MeshData<MSPACE>>::getInstance().at(id);
//   }

//   template <MemorySpace

// private:
//   MeshManager() = default;
//   ~MeshManager() = default;
// };

} // end namespace tribol

#endif /* SRC_MESH_MESHMANAGER_HPP_ */