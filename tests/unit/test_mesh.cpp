#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

#include "mesh.hpp"

using Catch::Matchers::WithinAbs;
using hydro::mesh::Index;
using hydro::mesh::kFacesPerCell;
using hydro::mesh::kHexFaceNodes;
using hydro::mesh::kNodesPerFace;
using hydro::mesh::kNoCell;
using hydro::mesh::Mesh;
using hydro::mesh::MeshMapping;
using hydro::mesh::StructuredMeshSpec;
using hydro::mesh::apply_saltzman_skew;
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
//   - geometrically, each face's right-hand normal points out of its owner
//     (skipped for a zero-area face, the x_min faces on a polar axis).
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
    if (normal == std::array<double, 3>{0.0, 0.0, 0.0}) continue;
    double outward = 0.0;
    for (int d = 0; d < 3; ++d) outward += normal[d] * (face_c[d] - cell_c[d]);
    CHECK(outward > 0.0);
  }
}

// Area of cell c's bottom (-z) quad by the shoelace formula; on a mesh
// extruded in z the cell's volume is this times its height.
double bottom_area(const Mesh& m, Index c) {
  const auto& cn = m.cell_nodes[c];
  double area = 0.0;
  for (int a = 0; a < 4; ++a) {
    const Index p = cn[a], q = cn[(a + 1) % 4];
    area += 0.5 * (m.node_x[p] * m.node_y[q] - m.node_x[q] * m.node_y[p]);
  }
  return area;
}

// Total area vector of face f (out of its owner).
std::array<double, 3> face_vector(const Mesh& m, Index f) {
  std::array<double, 3> total{};
  for (const auto& v : face_area_vectors(m, f)) {
    for (int d = 0; d < 3; ++d) total[d] += v[d];
  }
  return total;
}

// GCL: each cell's outward face vectors sum to zero, relative to its largest.
void check_closed_cells(const Mesh& m) {
  for (Index c = 0; c < m.num_cells(); ++c) {
    std::array<double, 3> sum{};
    double scale = 0.0;
    for (int l = 0; l < kFacesPerCell; ++l) {
      const Index f = m.cell_faces[c][l];
      const double sign = m.face_owner[f] == c ? 1.0 : -1.0;
      const auto v = face_vector(m, f);
      for (int d = 0; d < 3; ++d) {
        sum[d] += sign * v[d];
        scale = std::max(scale, std::abs(v[d]));
      }
    }
    INFO("cell " << c);
    for (int d = 0; d < 3; ++d) CHECK_THAT(sum[d], WithinAbs(0.0, 1e-15 * scale));
  }
}

StructuredMeshSpec polar_spec(Index nx, Index ny, Index nz, double r_min, double r_max) {
  StructuredMeshSpec spec{nx, ny, nz, r_min, r_max, 0.0, 90.0, 0.0, 0.1};
  spec.mapping = MeshMapping::Polar;
  return spec;
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

TEST_CASE("generate_structured_mesh: boundary faces tagged with their side",
          "[mesh]") {
  StructuredMeshSpec spec;
  spec.nx = 3;
  spec.ny = 2;
  spec.nz = 2;
  spec.x_max = 3.0;
  spec.y_max = 1.0;
  spec.z_max = 0.5;
  const Mesh m = generate_structured_mesh(spec);
  REQUIRE(m.face_boundary.size() == static_cast<std::size_t>(m.num_faces()));

  // Side s is local face s of the owner (kHexFaceNodes order), so a boundary
  // face's center lies on that side's plane.
  const std::array<double, kFacesPerCell> plane = {
      spec.x_min, spec.x_max, spec.y_min, spec.y_max, spec.z_min, spec.z_max};
  std::array<Index, kFacesPerCell> count{};
  for (Index f = 0; f < m.num_faces(); ++f) {
    INFO("face " << f);
    const int side = m.face_boundary[f];
    if (m.face_neighbor[f] != kNoCell) {
      CHECK(side == hydro::mesh::kInteriorFace);
      continue;
    }
    REQUIRE(side >= 0);
    REQUIRE(side < kFacesPerCell);
    ++count[side];
    CHECK(m.cell_faces[m.face_owner[f]][side] == f);
    double center = 0.0;
    for (Index n : m.face_nodes[f]) center += node(m, n)[side / 2] / kNodesPerFace;
    CHECK_THAT(center, WithinAbs(plane[side], 1e-14));
  }
  CHECK(count == std::array<Index, kFacesPerCell>{4, 4, 6, 6, 6, 6});
}

TEST_CASE("generate_structured_mesh: x_nodes grade the mesh in x", "[mesh]") {
  // 3x2x1 with x nodes 0, 0.1, 0.5, 1: same connectivity as the uniform
  // mesh, x from x_nodes at every (j, k), y and z unchanged.
  StructuredMeshSpec spec{3, 2, 1, 0.0, 1.0, 0.0, 1.0, 0.0, 0.5};
  const Mesh uniform = generate_structured_mesh(spec);
  spec.x_nodes = {0.0, 0.1, 0.5, 1.0};
  const Mesh m = generate_structured_mesh(spec);
  CHECK(m.cell_nodes == uniform.cell_nodes);
  CHECK(m.cell_faces == uniform.cell_faces);
  CHECK(m.face_nodes == uniform.face_nodes);
  CHECK(m.face_boundary == uniform.face_boundary);
  CHECK(m.node_y == uniform.node_y);
  CHECK(m.node_z == uniform.node_z);
  for (Index p = 0; p < m.num_nodes(); ++p) {
    CHECK(m.node_x[p] == spec.x_nodes[p % 4]);
  }
  check_face_consistency(m);
  spec.x_nodes = {0.0, 0.5, 1.0};
  CHECK_THROWS(generate_structured_mesh(spec));
}

TEST_CASE("apply_saltzman_skew: the standard 100x10 Saltzman mesh", "[mesh]") {
  // x_ij = i dx + (10 - j) dy sin(pi i dx), dx = dy = 0.01; z unchanged.
  StructuredMeshSpec spec{100, 10, 1, 0.0, 1.0, 0.0, 0.1, 0.0, 0.01};
  Mesh m = generate_structured_mesh(spec);
  const Mesh plain = m;
  apply_saltzman_skew(m, spec);
  auto node = [](Index i, Index j, Index k) { return k * 11 * 101 + j * 101 + i; };
  for (Index k = 0; k <= 1; ++k) {
    CHECK_THAT(m.node_x[node(50, 0, k)], WithinAbs(0.6, 1e-15));
    CHECK_THAT(m.node_x[node(50, 10, k)], WithinAbs(0.5, 1e-15));
    CHECK_THAT(m.node_x[node(25, 4, k)],
               WithinAbs(0.25 + 0.06 * std::sin(std::acos(-1.0) / 4.0), 1e-15));
    for (Index j = 0; j <= 10; ++j) {
      CHECK_THAT(m.node_x[node(0, j, k)], WithinAbs(0.0, 1e-15));
      CHECK_THAT(m.node_x[node(100, j, k)], WithinAbs(1.0, 1e-15));
    }
  }
  CHECK(m.node_y == plain.node_y);
  CHECK(m.node_z == plain.node_z);

  // Every cell stays valid: its bottom quad (planar in z) has positive area,
  // and the areas add up to the unchanged domain.
  double total = 0.0;
  for (Index c = 0; c < m.num_cells(); ++c) {
    const auto& cn = m.cell_nodes[c];
    double area = 0.0;
    for (int a = 0; a < 4; ++a) {
      const Index p = cn[a], q = cn[(a + 1) % 4];
      area += 0.5 * (m.node_x[p] * m.node_y[q] - m.node_x[q] * m.node_y[p]);
    }
    INFO("cell " << c);
    CHECK(area > 0.0);
    total += area;
  }
  CHECK_THAT(total, WithinAbs(0.1, 1e-14));
}

TEST_CASE("polar mapping: a quarter annulus", "[mesh]") {
  // r in [0.5, 1] (4 cells), theta in [0, 90] degrees (6 cells), z in
  // [0, 0.1] (2 cells): the structured numbering and connectivity, nodes at
  // (r cos theta, r sin theta, z).
  const StructuredMeshSpec spec = polar_spec(4, 6, 2, 0.5, 1.0);
  const Mesh m = generate_structured_mesh(spec);
  StructuredMeshSpec box = spec;
  box.mapping = MeshMapping::Cartesian;
  const Mesh cartesian = generate_structured_mesh(box);
  CHECK(m.cell_nodes == cartesian.cell_nodes);
  CHECK(m.face_nodes == cartesian.face_nodes);
  CHECK(m.face_boundary == cartesian.face_boundary);
  CHECK(m.node_z == cartesian.node_z);
  const double pi = std::acos(-1.0);
  for (Index p = 0; p < m.num_nodes(); ++p) {
    const double r = cartesian.node_x[p], theta = cartesian.node_y[p] * pi / 180.0;
    INFO("node " << p);
    CHECK_THAT(m.node_x[p], WithinAbs(r * std::cos(theta), 1e-15));
    CHECK_THAT(m.node_y[p], WithinAbs(r * std::sin(theta), 1e-15));
  }
  auto node = [](Index i, Index j, Index k) { return k * 7 * 5 + j * 5 + i; };
  for (Index k = 0; k <= 2; ++k) {
    for (Index i = 0; i <= 4; ++i) {
      // The planes theta = 0 and 90 exactly; theta = 30 and 60 mirror
      // images about the diagonal, bit for bit.
      CHECK(m.node_y[node(i, 0, k)] == 0.0);
      CHECK(m.node_x[node(i, 6, k)] == 0.0);
      CHECK(m.node_x[node(i, 2, k)] == m.node_y[node(i, 4, k)]);
      CHECK(m.node_y[node(i, 2, k)] == m.node_x[node(i, 4, k)]);
    }
  }
  check_face_consistency(m);
  check_closed_cells(m);

  // Every cell is an extruded planar quad, 1/2 (r1^2 - r0^2) sin(dtheta)
  // in area; the total is the polygonal approximation of the sector.
  const double dtheta = pi / 12.0;
  double total = 0.0;
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Index i = c % 4;
    const double r0 = 0.5 + 0.125 * i, r1 = r0 + 0.125;
    INFO("cell " << c);
    CHECK_THAT(bottom_area(m, c), WithinAbs(0.5 * (r1 * r1 - r0 * r0) * std::sin(dtheta), 1e-15));
    total += bottom_area(m, c) * 0.05;
  }
  const double polygon = 6 * 0.5 * (1.0 - 0.25) * std::sin(dtheta) * 0.1;
  CHECK_THAT(total, WithinAbs(polygon, 1e-15));
  // ... which approaches the annular sector's pi/4 (1 - 0.25) 0.1 as
  // sin(dtheta)/dtheta -> 1.
  CHECK_THAT(total / (pi / 4.0 * 0.75 * 0.1), WithinAbs(std::sin(dtheta) / dtheta, 1e-14));
}

TEST_CASE("polar mapping: r_min = 0 has one node per level on the axis", "[mesh]") {
  // 3x4x2 down to r = 0: the axis nodes (0, j, k) are one node per k, the
  // i = 0 cells triangular prisms, their x_min faces of zero area.
  const Mesh m = generate_structured_mesh(polar_spec(3, 4, 2, 0.0, 1.0));
  CHECK(m.num_nodes() == (3 * 5 + 1) * 3);
  CHECK(m.num_cells() == 3 * 4 * 2);
  // Nodes: per level, the axis node, then (i >= 1, j = 0), then (i >= 1, j)
  // for j = 1..4; level k starts at 16 k.
  auto node = [](Index i, Index j, Index k) {
    return i == 0 ? 16 * k : 16 * k + 1 + 3 * j + (i - 1);
  };
  for (Index k = 0; k <= 2; ++k) {
    const Index a = node(0, 0, k);
    CHECK(m.node_x[a] == 0.0);
    CHECK(m.node_y[a] == 0.0);
    CHECK_THAT(m.node_z[a], WithinAbs(0.05 * k, 1e-15));
    CHECK(m.node_x[node(2, 4, k)] == 0.0);
    CHECK_THAT(m.node_y[node(2, 4, k)], WithinAbs(2.0 / 3.0, 1e-15));
  }
  for (Index k = 0; k < 2; ++k) {
    for (Index j = 0; j < 4; ++j) {
      const auto& cn = m.cell_nodes[k * 12 + j * 3];
      INFO("axis cell j = " << j << ", k = " << k);
      CHECK(cn[0] == node(0, 0, k));
      CHECK(cn[3] == node(0, 0, k));
      CHECK(cn[4] == node(0, 0, k + 1));
      CHECK(cn[7] == node(0, 0, k + 1));
      CHECK(cn[1] == node(1, j, k));
      CHECK(cn[2] == node(1, j + 1, k));
    }
  }
  // Every node is in some cell.
  std::vector<int> uses(m.num_nodes(), 0);
  for (const auto& cn : m.cell_nodes) {
    for (const Index p : cn) ++uses[p];
  }
  CHECK(std::count(uses.begin(), uses.end(), 0) == 0);

  // Zero area exactly on the axis faces, nowhere else.
  for (Index f = 0; f < m.num_faces(); ++f) {
    const auto v = face_vector(m, f);
    const bool zero = v == std::array<double, 3>{0.0, 0.0, 0.0};
    INFO("face " << f);
    CHECK(zero == (m.face_boundary[f] == 0));
  }
  check_face_consistency(m);
  check_closed_cells(m);
  const double pi = std::acos(-1.0);
  for (Index c = 0; c < m.num_cells(); ++c) {
    const Index i = c % 3;
    const double r0 = i / 3.0, r1 = (i + 1) / 3.0;
    INFO("cell " << c);
    CHECK_THAT(bottom_area(m, c),
               WithinAbs(0.5 * (r1 * r1 - r0 * r0) * std::sin(pi / 8.0), 1e-15));
  }
}

TEST_CASE("polar mapping: rejected ranges", "[mesh]") {
  CHECK_THROWS(generate_structured_mesh(polar_spec(2, 2, 1, -0.1, 1.0)));
  StructuredMeshSpec spec = polar_spec(2, 2, 1, 0.5, 1.0);
  spec.y_max = 180.5;
  CHECK_THROWS(generate_structured_mesh(spec));
  spec.y_min = spec.y_max = 10.0;
  CHECK_THROWS(generate_structured_mesh(spec));
  // A half plane is fine in two cells, not in one.
  spec.y_min = -90.0;
  spec.y_max = 90.0;
  const Mesh half = generate_structured_mesh(spec);
  CHECK(half.node_x[0] == 0.0);
  CHECK(half.node_y[0] == -0.5);
  check_closed_cells(half);
  spec.ny = 1;
  CHECK_THROWS(generate_structured_mesh(spec));
}
