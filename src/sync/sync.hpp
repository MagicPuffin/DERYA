#pragma once

namespace hydro::sync {

// The MPI-aware layer (decisions.md, "Parallelism"). Until v0.4 there is a
// single rank and every operation is the identity.
enum class ReduceOp { Min, Max, Sum };

// Reduction of one value over all ranks.
inline double global_reduce(double local, ReduceOp /*op*/) { return local; }

}  // namespace hydro::sync
