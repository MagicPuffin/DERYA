#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <vector>

#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using hydro::mesh::Index;
using hydro::mesh::kFacesPerCell;
using hydro::mesh::kHexFaceNodes;
using hydro::mesh::kNodesPerFace;
using hydro::mesh::kNoCell;
using hydro::mesh::Mesh;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::face_area_vectors;
using hydro::mesh::generate_structured_mesh;

namespace {

using Quad = std::array<Index, kNodesPerFace>;

// True if b is a cyclic rotation of a.
bool same_cycle(const Quad& a, const Quad& b) {
  for (int s = 0; s < kNodesPerFace; ++s) {
    bool match = true;
    for (int n = 0; n < kNodesPerFace; ++n) {
      if (a[(n + s) % kNodesPerFace] != b[n]) match = false;
    }
    if (match) return true;
  }
  return false;
}

std::array<double, 3> node(const Mesh& m, Index n) {
  return {m.node_x[n], m.node_y[n], m.node_z[n]};
}

// Structural invariants that must hold for any mesh, independent of numbering:
//   - each cell's local face l, read through kHexFaceNodes, is face
//     cell_faces[l] with the same orientation if the cell owns it and the
//     reversed orientation if it is the neighbor;
//   - geometrically, each face's right-hand normal points out of its owner.
void check_face_consistency(const Mesh& m) {
  for (Index c = 0; c < m.num_cells(); ++c) {
    for (int l = 0; l < kFacesPerCell; ++l) {
      INFO("cell " << c << ", local face " << l);
      const Index f = m.cell_faces[c][l];
      Quad local;
      for (int n = 0; n < kNodesPerFace; ++n) {
        local[n] = m.cell_nodes[c][kHexFaceNodes[l][n]];
      }
      if (m.face_owner[f] == c) {
        CHECK(same_cycle(local, m.face_nodes[f]));
      } else {
        CHECK(m.face_neighbor[f] == c);
        const Quad reversed = {local[3], local[2], local[1], local[0]};
        CHECK(same_cycle(reversed, m.face_nodes[f]));
      }
    }
  }

  for (Index f = 0; f < m.num_faces(); ++f) {
    INFO("face " << f);
    const auto& fn = m.face_nodes[f];
    const auto p0 = node(m, fn[0]), p1 = node(m, fn[1]);
    const auto p2 = node(m, fn[2]), p3 = node(m, fn[3]);
    // Normal from the diagonals: (p2 - p0) x (p3 - p1).
    const std::array<double, 3> d1 = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
    const std::array<double, 3> d2 = {p3[0] - p1[0], p3[1] - p1[1], p3[2] - p1[2]};
    const std::array<double, 3> normal = {d1[1] * d2[2] - d1[2] * d2[1],
                                          d1[2] * d2[0] - d1[0] * d2[2],
                                          d1[0] * d2[1] - d1[1] * d2[0]};
    std::array<double, 3> face_c{}, cell_c{};
    for (int n = 0; n < kNodesPerFace; ++n) {
      const auto p = node(m, fn[n]);
      for (int d = 0; d < 3; ++d) face_c[d] += p[d] / kNodesPerFace;
    }
    for (Index n : m.cell_nodes[m.face_owner[f]]) {
      const auto p = node(m, n);
      for (int d = 0; d < 3; ++d) cell_c[d] += p[d] / 8.0;
    }
    double outward = 0.0;
    for (int d = 0; d < 3; ++d) outward += normal[d] * (face_c[d] - cell_c[d]);
    CHECK(outward > 0.0);
  }
}

}  // namespace

TEST_CASE("generate_structured_mesh: 2x2x1 slab counts and cell 0 corners",
          "[mesh]") {
  StructuredMeshSpec spec;
  spec.nx = 2;
  spec.ny = 2;  // nz left at its default of 1
  const Mesh m = generate_structured_mesh(spec);

  CHECK(m.num_nodes() == 18);
  CHECK(m.num_cells() == 4);
  CHECK(m.num_faces() == 20);
  REQUIRE(m.num_cells() >= 1);
  CHECK(m.cell_nodes[0] == std::array<Index, 8>{0, 1, 4, 3, 9, 10, 13, 12});

  REQUIRE(m.node_y.size() == 18);
  REQUIRE(m.node_z.size() == 18);
  REQUIRE(m.cell_faces.size() == 4);
  REQUIRE(m.face_owner.size() == 20);
  REQUIRE(m.face_neighbor.size() == 20);
  check_face_consistency(m);
}

TEST_CASE("generate_structured_mesh: 3x2x2 counts and consistency", "[mesh]") {
  // Distinct nx/ny/nz and a non-cubic domain, so swapped axes show up.
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.x_max = 3.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  const Mesh m = generate_structured_mesh(spec);

  REQUIRE(m.num_nodes() == 4 * 3 * 3);
  REQUIRE(m.num_cells() == 3 * 2 * 2);
  REQUIRE(m.num_faces() == 4 * 2 * 2 + 3 * 3 * 2 + 3 * 2 * 3);

  Index boundary_faces = 0;
  for (Index f = 0; f < m.num_faces(); ++f) {
    if (m.face_neighbor[f] == kNoCell) ++boundary_faces;
  }
  CHECK(boundary_faces == 2 * (2 * 2 + 3 * 2 + 3 * 2));

  check_face_consistency(m);
}

namespace {

using Vec = std::array<double, 3>;

void check_vec(const Vec& actual, const Vec& expected) {
  for (int d = 0; d < 3; ++d) CHECK_THAT(actual[d], WithinAbs(expected[d], 1e-14));
}

// Sum over a cell's faces and nodes of the outward S_pf n_pf, and the cell
// volume from them via V_c = 1/3 sum_p x_p . sum_f S_pf n_pf (exact for the
// face-split cell, since V_c is cubic in the node positions and the GCL gives
// dV_c/dx_p = sum_f S_pf n_pf).
void cell_closure_and_volume(const Mesh& m, Index c, Vec& closure, double& volume) {
  closure = {0.0, 0.0, 0.0};
  volume = 0.0;
  for (Index f : m.cell_faces[c]) {
    const double sign = m.face_owner[f] == c ? 1.0 : -1.0;
    const auto vecs = face_area_vectors(m, f);
    for (int n = 0; n < kNodesPerFace; ++n) {
      const auto x = node(m, m.face_nodes[f][n]);
      for (int d = 0; d < 3; ++d) {
        closure[d] += sign * vecs[n][d];
        volume += sign * x[d] * vecs[n][d] / 3.0;
      }
    }
  }
}

}  // namespace

TEST_CASE("face_area_vectors: axis-aligned faces of a generated cell", "[mesh]") {
  // dx = 1, dy = 0.5, dz = 0.25: every face of a rectangular cell splits
  // evenly, so each node gets a quarter of the face area along the normal.
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.x_max = 3.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  const Mesh m = generate_structured_mesh(spec);

  // Cell 0 owns all six of its faces (boundary, or lower-indexed neighbor).
  const Index c = 0;
  const std::array<Vec, 6> expected = {{
      {-0.125 / 4, 0.0, 0.0},  // -x: area dy*dz = 0.125
      {+0.125 / 4, 0.0, 0.0},  // +x
      {0.0, -0.25 / 4, 0.0},   // -y: area dx*dz = 0.25
      {0.0, +0.25 / 4, 0.0},   // +y
      {0.0, 0.0, -0.5 / 4},    // -z: area dx*dy = 0.5
      {0.0, 0.0, +0.5 / 4},    // +z
  }};
  for (int l = 0; l < kFacesPerCell; ++l) {
    const Index f = m.cell_faces[c][l];
    REQUIRE(m.face_owner[f] == c);
    const auto vecs = face_area_vectors(m, f);
    for (int n = 0; n < kNodesPerFace; ++n) {
      INFO("local face " << l << ", face node " << n);
      check_vec(vecs[n], expected[l]);
    }
  }

  Vec closure;
  double volume;
  cell_closure_and_volume(m, c, closure, volume);
  check_vec(closure, {0.0, 0.0, 0.0});
  CHECK_THAT(volume, WithinAbs(1.0 * 0.5 * 0.25, 1e-14));
}

TEST_CASE("face_area_vectors: trapezoidal face weights nodes unequally", "[mesh]") {
  // Unit-height cell on [0,2]x[0,1]x[0,1] with bottom node 3 moved from
  // (2,1,0) to (1,1,0): the -z face becomes the trapezoid (0,0), (2,0), (1,1),
  // (0,1), area 1.5, barycenter p* = (0.75, 0.5).
  StructuredMeshSpec spec;
  spec.x_max = 2.0;
  Mesh m = generate_structured_mesh(spec);
  m.node_x[3] = 1.0;

  // Triangles (p*, p, p+) around the trapezoid, areas by hand:
  //   (0,0)-(2,0): 0.5    (2,0)-(1,1): 0.375
  //   (1,1)-(0,1): 0.25   (0,1)-(0,0): 0.375
  // Eq. 5, node p: 1/3 (its two triangles + 1.5/4)
  //   (0,0), (2,0): 1/3 (0.875 + 0.375) = 5/12
  //   (1,1), (0,1): 1/3 (0.625 + 0.375) = 1/3
  // (An even split would give 0.375 each.) The face normal is -z.
  const Index f = m.cell_faces[0][4];  // local face -z, nodes {0, 2, 3, 1}
  REQUIRE(m.face_nodes[f] == Quad{0, 2, 3, 1});
  const auto vecs = face_area_vectors(m, f);
  check_vec(vecs[0], {0.0, 0.0, -5.0 / 12.0});  // node 0 (0,0,0)
  check_vec(vecs[1], {0.0, 0.0, -1.0 / 3.0});   // node 2 (0,1,0)
  check_vec(vecs[2], {0.0, 0.0, -1.0 / 3.0});   // node 3 (1,1,0)
  check_vec(vecs[3], {0.0, 0.0, -5.0 / 12.0});  // node 1 (2,0,0)

  // The +x and +y faces are now non-planar. Volume by the prismatoid formula,
  // h/6 (A_bottom + 4 A_mid + A_top), mid section a trapezoid with top edge
  // 1.5: 1/6 (1.5 + 4*1.75 + 2) = 1.75.
  Vec closure;
  double volume;
  cell_closure_and_volume(m, 0, closure, volume);
  check_vec(closure, {0.0, 0.0, 0.0});
  CHECK_THAT(volume, WithinAbs(1.75, 1e-14));
}
