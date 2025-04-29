#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <vector>
#include <map>
#include <string>
#include <tuple>
#include <algorithm>
#include "utils.h"  // for utils::vec and utils::PBCMask

namespace geometry {

using vec = utils::vec;

// ===== Structures =====

struct am_interaction {
    double myosin_binding_ratio = 0.0;
    vec myosin_binding_start;
    vec myosin_binding_end;
    vec crosslinkable_start;
    vec crosslinkable_end;
    double partial_binding_ratio = 0.0;
    double crosslinkable_ratio = 0.0;
};

// ===== Core Functions =====

// Compute the shortest distance between a point and a segment (no PBC)
double point_segment_distance(const vec& x, const vec& a, const vec& b);

// Compute the shortest distance and normal vector between a point and a segment (no PBC)
std::pair<double, vec> point_segment_distance_w_normal(const vec& x, const vec& a, const vec& b);

// Compute the orientation of three points (no PBC)
int orientation(const vec& p, const vec& q, const vec& r);

// Apply PBC between four points (shift actin to match myosin) — modifies arguments
void apply_pbc(
    vec& actin_left, vec& actin_right,
    vec& myosin_left, vec& myosin_right,
    const std::vector<double>& box,
    const utils::PBCMask& pbc_mask
);

// Compute the shortest distance between two segments (with PBC)
double segment_segment_distance(
    const vec& a, const vec& b, const vec& c, const vec& d,
    const std::vector<double>& box,
    const utils::PBCMask& pbc_mask
);

// Compute the shortest distance between two segments (with PBC)
double segment_segment_distance(
    const vec& a, const vec& b, const vec& c, const vec& d
);
// Compute the shortest distance and normal vector between two segments (with PBC)
std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(
    const vec& a, const vec& b, const vec& c, const vec& d,
    const std::vector<double>& box,
    const utils::PBCMask& pbc_mask
);

// ===== Subsegment Finding =====

// Solve a quadratic equation: a * t^2 + b * t + c = 0
std::vector<double> solveQuadratic(double a, double b, double c);

// Find roots of a quadratic that lie within [t_start, t_end]
std::vector<double> findRootsInInterval(const std::vector<double>& roots, double t_start, double t_end);

// Find the subsegment on AB where distance to CD is < d
std::tuple<double, vec, vec> subsegment_within_distance(
    const vec& A, const vec& B,
    const vec& C, const vec& D,
    double& d
);

// Same as above, but applies PBC before measuring
std::tuple<double, vec, vec> subsegment_within_distance(
    const vec& A, const vec& B, const vec& C, const vec& D,
    double& d,
    const std::vector<double>& box,
    const utils::PBCMask& pbc_mask
);

// ===== Main Interaction Analysis =====

// Analyze actin–myosin geometry and return interaction details
am_interaction analyze_am(
    const vec& actin_left, const vec& actin_right,
    const vec& myosin_left, const vec& myosin_right,
    double& d,
    const std::vector<double>& box,
    const utils::PBCMask& pbc_mask
);

} // namespace geometry

#endif  // GEOMETRY_H
