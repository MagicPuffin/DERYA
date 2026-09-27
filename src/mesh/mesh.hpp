#pragma once

#include <array>
#include <cstdint>
#include <vector>

namespace hydro::mesh {

// Local entity index. Plain integer for now; the local-index-plus-ghost-flag
// type arrives with the Task 4 index types.
using Index = std::int32_t;

// Face neighbor of a boundary face.
inline constexpr Index kNoCell = -1;

inline constexpr int kNodesPerCell = 8;
inline constexpr int kFacesPerCell = 6;
inline constexpr int kNodesPerFace = 4;

// Local face l of a hex, as local corner indices (positions in cell_nodes).
// Faces are ordered -x, +x, -y, +y, -z, +z; each face's nodes follow the
// right-hand rule with the normal pointing out of the cell.
inline constexpr std::array<std::array<int, kNodesPerFace>, kFacesPerCell>
    kHexFaceNodes = {{
        {0, 4, 7, 3},  // -x
        {1, 2, 6, 5},  // +x
        {0, 1, 5, 4},  // -y
        {3, 7, 6, 2},  // +y
        {0, 3, 2, 1},  // -z
        {4, 5, 6, 7},  // +z
    }};

// 3D unstructured hexahedral mesh: node coordinates (SoA) plus connectivity.
// "2D" problems run as a one-cell-thick slab (nz = 1), see decisions.md.
//
// Conventions (the structured generator and everything downstream rely on
// these):
//   - cell_nodes lists a hex's 8 corners: the bottom (-z) quad
//     counter-clockwise seen from +z, then the top quad in the same order.
//   - cell_faces[l] is the face on local face l of kHexFaceNodes.
//   - Every face has an owner cell; face_nodes is ordered so the right-hand
//     normal points out of the owner, into the neighbor. Boundary faces have
//     neighbor kNoCell.
struct Mesh {
  std::vector<double> node_x;
  std::vector<double> node_y;
  std::vector<double> node_z;

  std::vector<std::array<Index, kNodesPerCell>> cell_nodes;
  std::vector<std::array<Index, kFacesPerCell>> cell_faces;

  std::vector<std::array<Index, kNodesPerFace>> face_nodes;
  std::vector<Index> face_owner;
  std::vector<Index> face_neighbor;

  Index num_nodes() const { return static_cast<Index>(node_x.size()); }
  Index num_cells() const { return static_cast<Index>(cell_nodes.size()); }
  Index num_faces() const { return static_cast<Index>(face_nodes.size()); }
};

// Uniform nx-by-ny-by-nz hex mesh on [x_min, x_max] x [y_min, y_max] x
// [z_min, z_max]. nz = 1 is the thin-slab "2D" mode.
struct StructuredMeshSpec {
  Index nx = 1;
  Index ny = 1;
  Index nz = 1;
  double x_min = 0.0;
  double x_max = 1.0;
  double y_min = 0.0;
  double y_max = 1.0;
  double z_min = 0.0;
  double z_max = 1.0;
};

// Numbering, with node (i, j, k) at (x_min + i*dx, y_min + j*dy, z_min + k*dz):
//   - node (i, j, k), 0 <= i <= nx, 0 <= j <= ny, 0 <= k <= nz:
//     k*(ny+1)*(nx+1) + j*(nx+1) + i
//   - cell (i, j, k), 0 <= i < nx, 0 <= j < ny, 0 <= k < nz:
//     k*ny*nx + j*nx + i, corners (i,j,k), (i+1,j,k), (i+1,j+1,k), (i,j+1,k),
//     then the same four at k+1
//   - x-normal faces first, face (i, j, k) on x = x_i, i <= nx:
//     k*ny*(nx+1) + j*(nx+1) + i
//   - then y-normal faces, face (i, j, k) on y = y_j, j <= ny:
//     Fx + k*(ny+1)*nx + j*nx + i,          Fx = (nx+1)*ny*nz
//   - then z-normal faces, face (i, j, k) on z = z_k, k <= nz:
//     Fx + Fy + k*ny*nx + j*nx + i,         Fy = nx*(ny+1)*nz
//   - an interior face is owned by the lower-indexed of its two cells.
Mesh generate_structured_mesh(const StructuredMeshSpec& spec);

// Nodal face area vectors S_pf n_pf of face f (Colaïtis, Guisset & Breil,
// SSRN 5167115, Eq. 5, after Georges, Breil & Maire [29]). The face is split
// into triangles tr = (p*_f, p, p+) around its barycenter p*_f, and node p gets
//   S_pf n_pf = 1/3 [ sum_{tr in T(f,p)} S_tr n_tr
//                     + 1/N_f sum_{tr in T(f)} S_tr n_tr ],
// with N_f the number of face nodes; the entries sum to the face area vector.
// Entry n belongs to node face_nodes[f][n]. Vectors point out of the owner
// cell; the neighbor cell uses their negatives.
std::array<std::array<double, 3>, kNodesPerFace> face_area_vectors(
    const Mesh& m, Index f);

}  // namespace hydro::mesh
