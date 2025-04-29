#include "geometry.h"
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <numeric>
#include <utility>
#include <tuple>
#include <string>

namespace geometry {

const double EPS = 1e-9;

// Helper: copy 4 vectors
void copy_segments(const vec& a, const vec& b, const vec& c, const vec& d, vec& a_copy, vec& b_copy, vec& c_copy, vec& d_copy) {
    a_copy = a;
    b_copy = b;
    c_copy = c;
    d_copy = d;
}

// ------------------ Distance and Orientation Functions ------------------

double point_segment_distance(const vec& x, const vec& a, const vec& b) {
    vec ab = b - a;
    double ab_norm = ab.norm();
    vec ab_normalized = {ab.x / ab_norm, ab.y / ab_norm};
    vec ap = x - a;
    double ap_ab = ap.x * ab_normalized.x + ap.y * ab_normalized.y;
    ap_ab = std::fmax(0, std::fmin(ap_ab, ab_norm));
    vec closest_point = a + vec{ap_ab * ab_normalized.x, ap_ab * ab_normalized.y};
    vec norm_vec = x - closest_point;
    return norm_vec.norm();
}

std::pair<double, vec> point_segment_distance_w_normal(const vec& x, const vec& a, const vec& b) {
    vec ab = b - a;
    double ab_norm = ab.norm();
    vec ab_normalized = {ab.x / ab_norm, ab.y / ab_norm};
    vec ap = x - a;
    double ap_ab = ap.x * ab_normalized.x + ap.y * ab_normalized.y;
    ap_ab = std::fmax(0, std::fmin(ap_ab, ab_norm));
    vec norm_vec = ap - vec{ap_ab * ab_normalized.x, ap_ab * ab_normalized.y};
    double dist = norm_vec.norm();
    return {dist, norm_vec};
}

int orientation(const vec& p, const vec& q, const vec& r) {
    vec pq = q - p;
    vec qr = r - q;
    double val = pq.y * qr.x - pq.x * qr.y;
    if (val == 0) return 0;
    return (val > 0) ? 1 : 2;
}

// ------------------ Apply PBC ------------------

void apply_pbc(vec& actin_left, vec& actin_right, vec& myosin_left, vec& myosin_right, const std::vector<double>& box, const utils::PBCMask& pbc_mask){
    vec actin_mid = (actin_left + actin_right) / 2;
    vec myosin_mid = (myosin_left + myosin_right) / 2;
    vec displacement = actin_mid - myosin_mid;
    vec shift;

    if (pbc_mask.x) shift.x = -box[0] * round(displacement.x / box[0]);
    else shift.x = 0;
    if (pbc_mask.y) shift.y = -box[1] * round(displacement.y / box[1]);
    else shift.y = 0;

    actin_left = actin_left + shift;
    actin_right = actin_right + shift;
}

// ------------------ Segment Distance Functions ------------------

double segment_segment_distance(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box, const utils::PBCMask& pbc_mask) {
    vec a_copy, b_copy, c_copy, d_copy;
    copy_segments(a, b, c, d, a_copy, b_copy, c_copy, d_copy);
    apply_pbc(a_copy, b_copy, c_copy, d_copy, box, pbc_mask);
    return segment_segment_distance(a_copy, b_copy, c_copy, d_copy);
}

double segment_segment_distance(const vec& a, const vec& b, const vec& c, const vec& d) {
    double a_cd = point_segment_distance(a, c, d);
    double b_cd = point_segment_distance(b, c, d);
    double c_ab = point_segment_distance(c, a, b);
    double d_ab = point_segment_distance(d, a, b);

    double dist = std::min({a_cd, b_cd, c_ab, d_ab});

    int o1 = orientation(a, b, c);
    int o2 = orientation(a, b, d);
    int o3 = orientation(c, d, a);
    int o4 = orientation(c, d, b);

    if ((o1 != o2) && (o3 != o4)) dist = 0;

    return dist;
}

std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box, const utils::PBCMask& pbc_mask) {
    vec a_copy, b_copy, c_copy, d_copy;
    copy_segments(a, b, c, d, a_copy, b_copy, c_copy, d_copy);
    apply_pbc(a_copy, b_copy, c_copy, d_copy, box, pbc_mask);

    auto a_cd = point_segment_distance_w_normal(a_copy, c_copy, d_copy);
    auto b_cd = point_segment_distance_w_normal(b_copy, c_copy, d_copy);
    auto c_ab = point_segment_distance_w_normal(c_copy, a_copy, b_copy);
    auto d_ab = point_segment_distance_w_normal(d_copy, a_copy, b_copy);

    double dist = std::min({a_cd.first, b_cd.first, c_ab.first, d_ab.first});
    std::map<std::string, vec> normal_dict;

    if (dist == a_cd.first) {
        normal_dict["start"] = a_copy;
        normal_dict["vector"] = a_cd.second;
        normal_dict["end"] = a_copy + a_cd.second;
    } else if (dist == b_cd.first) {
        normal_dict["start"] = b_copy;
        normal_dict["vector"] = b_cd.second;
        normal_dict["end"] = b_copy + b_cd.second;
    } else if (dist == c_ab.first) {
        normal_dict["end"] = c_copy;
        normal_dict["vector"] = -c_ab.second;
        normal_dict["start"] = c_copy - c_ab.second;
    } else {
        normal_dict["end"] = d_copy;
        normal_dict["vector"] = -d_ab.second;
        normal_dict["start"] = d_copy - d_ab.second;
    }

    int o1 = orientation(a_copy, b_copy, c_copy);
    int o2 = orientation(a_copy, b_copy, d_copy);
    int o3 = orientation(c_copy, d_copy, a_copy);
    int o4 = orientation(c_copy, d_copy, b_copy);

    if ((o1 != o2) && (o3 != o4)) {
        dist = 0;
        normal_dict["vector"] = {0, 0};
    }

    return {dist, normal_dict};
}

// ------------------ Root Finding ------------------

std::vector<double> findRootsInInterval(const std::vector<double>& roots, double t_start, double t_end) {
    std::vector<double> valid_roots;
    for (double t : roots) {
        if (t >= t_start - EPS && t <= t_end + EPS) valid_roots.push_back(t);
    }
    return valid_roots;
}

std::vector<double> solveQuadratic(double a, double b, double c) {
    std::vector<double> roots;
    roots.reserve(2);
    double discriminant = b * b - 4 * a * c;

    if (discriminant < -EPS) return roots;
    else if (std::abs(discriminant) < EPS) {
        roots.push_back(-b / (2 * a));
    } else {
        double sqrt_disc = std::sqrt(discriminant);
        roots.push_back((-b + sqrt_disc) / (2 * a));
        roots.push_back((-b - sqrt_disc) / (2 * a));
    }
    return roots;
}

// ------------------ Subsegment Search ------------------

std::tuple<double, vec, vec> subsegment_within_distance(const vec& A, const vec& B, const vec& C, const vec& D, double& d) {
    double min_dist = segment_segment_distance(A, B, C, D);
    if (min_dist > d) {
        return std::make_tuple(0.0, A, A);
    }
    std::vector<std::pair<vec, vec>> intervals;
    double interval_start = 1.0;
    double interval_end = 0.0;
    vec dirAB = B - A; // Direction vector of AB
    vec dirCD = D - C; // Direction vector of CD
    double lenAB_sq = dirAB.norm_squared();
    double lenCD_sq = dirCD.norm_squared();

    //rule out the case where the midpoints of the two segments are more than d apart
    vec AB_mid = (A + B) / 2;
    vec CD_mid = (C + D) / 2;
    if ((AB_mid - CD_mid).norm() > dirAB.norm() + dirCD.norm() + d) {
        return std::make_tuple(0.0, A, A);
    }
    // Find transition points where closest point on CD changes
    std::vector<double> boundaries = {0.0, 1.0};
    boundaries.reserve(4);
    if (lenAB_sq > EPS) {
        // Compute t where normals from C and D intersect AB
        vec AC = A - C;
        vec AD = A - D;

        double t_c = - (dirAB.dot(AC)) / lenAB_sq;
        double t_d = - (dirAB.dot(AD)) / lenAB_sq;

        if (t_c > EPS && t_c < 1.0 - EPS)
            boundaries.push_back(t_c);
        if (t_d > EPS && t_d < 1.0 - EPS)
            boundaries.push_back(t_d);
    }

    // Remove duplicates and sort boundaries
    std::sort(boundaries.begin(), boundaries.end());
    boundaries.erase(std::unique(boundaries.begin(), boundaries.end(), [](double a, double b) {
        return std::abs(a - b) < EPS;
    }), boundaries.end());

    // Process each subsegment
    for (size_t i = 0; i < boundaries.size() - 1; ++i) {
        double t_start = boundaries[i];
        double t_end = boundaries[i + 1];

        if (t_end - t_start < EPS) {
            continue; // Skip negligible intervals
        }

        // Determine the closest point on CD for this subsegment
        // Use midpoint to determine whether closest point is along CD or at an endpoint
        double t_mid = (t_start + t_end) / 2.0;
        vec P_mid = A + dirAB * t_mid;

        // Compute s for the projection of P_mid onto CD
        double s = (P_mid - C).dot(dirCD) / lenCD_sq;

        vec closest_point;
        bool closest_is_projection = true;

        if (s <= 0.0) {
            closest_point = C;
            closest_is_projection = false;
        } else if (s >= 1.0) {
            closest_point = D;
            closest_is_projection = false;
        } else {
            closest_point = C + dirCD * s;
        }

        if (closest_is_projection) {
            // Closest point lies along CD
            // Set up the quadratic equation for D(t)^2 = d^2
            vec E = A - C;
            double A_coef = lenCD_sq * dirAB.dot(dirAB) - pow(dirAB.dot(dirCD), 2);
            double B_coef = 2 * (lenCD_sq * dirAB.dot(E) - dirAB.dot(dirCD) * dirCD.dot(E));
            double C_coef = lenCD_sq * E.dot(E) - pow(dirCD.dot(E), 2) - lenCD_sq * d * d;

            // Solve quadratic equation A_coef * t^2 + B_coef * t + C_coef = 0
            std::vector<double> roots = solveQuadratic(A_coef, B_coef, C_coef);

            // Filter roots within [t_start, t_end]
            std::vector<double> valid_roots;
            valid_roots.reserve(roots.size());
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
            interval_points.reserve(valid_roots.size() + 2);
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);

            // Sort and remove duplicates
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            // Determine intervals where D(t) < d
            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];

                if (t1 - t0 < EPS) {
                    continue;
                }

                // Test midpoint
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;

                // Compute s for P_test
                double s_test = (P_test - C).dot(dirCD) / lenCD_sq;
                s_test = std::max(0.0, std::min(1.0, s_test));
                vec closest_test = C + dirCD * s_test;
                double dist_sq = (P_test - closest_test).norm_squared();

                if (dist_sq < d * d + EPS) {
                    // Interval where distance < d
                    vec start_point = A + dirAB * t0;
                    vec end_point = A + dirAB * t1;
                    if (t0<interval_start){
                        interval_start = t0;
                    }
                    if (t1>interval_end){
                        interval_end = t1;
                    }
                    intervals.push_back({start_point, end_point});
                }
            }
        } else {
            // Closest point is at C or D
            vec Q = closest_point;
            vec M = A - Q;
            vec N = dirAB;

            double a_quad = N.dot(N);
            double b_quad = 2 * M.dot(N);
            double c_quad = M.dot(M) - d * d;

            // Solve quadratic equation a_quad * t^2 + b_quad * t + c_quad = 0
            std::vector<double> roots = solveQuadratic(a_quad, b_quad, c_quad);

            // Filter roots within [t_start, t_end]
            std::vector<double> valid_roots;
            valid_roots.reserve(roots.size());
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
            interval_points.reserve(valid_roots.size() + 2);
            interval_points.insert(interval_points.end(), valid_roots.begin(), valid_roots.end());
            interval_points.push_back(t_end);

            // Sort and remove duplicates
            std::sort(interval_points.begin(), interval_points.end());
            interval_points.erase(std::unique(interval_points.begin(), interval_points.end(), [](double a, double b) {
                return std::abs(a - b) < EPS;
            }), interval_points.end());

            // Determine intervals where D(t) < d
            for (size_t j = 0; j < interval_points.size() - 1; ++j) {
                double t0 = interval_points[j];
                double t1 = interval_points[j + 1];

                if (t1 - t0 < EPS) {
                    continue;
                }

                // Test midpoint
                double t_test = (t0 + t1) / 2.0;
                vec P_test = A + dirAB * t_test;
                double dist_sq = (P_test - Q).norm_squared();

                if (dist_sq < d * d + EPS) {
                    // Interval where distance < d
                    vec start_point = A + dirAB * t0;
                    vec end_point = A + dirAB * t1;
                    intervals.push_back({start_point, end_point});
                    if (t0<interval_start){
                        interval_start = t0;
                    }
                    if (t1>interval_end){
                        interval_end = t1;
                    }
                }
            }
        }
    }
    if (interval_start<=interval_end){
        vec inter = A + dirAB * interval_start;
        return {interval_end-interval_start, A + dirAB * interval_start, A+dirAB*interval_end};
    }
    else{
        return {0.0,A,A};
    }
}


// Wrapper with PBC
std::tuple<double, vec, vec> subsegment_within_distance(const vec& A, const vec& B, const vec& C, const vec& D, double& d, const std::vector<double>& box, const utils::PBCMask& pbc_mask) {
    vec A_copy = A, B_copy = B, C_copy = C, D_copy = D;
    apply_pbc(A_copy, B_copy, C_copy, D_copy, box, pbc_mask);
    return subsegment_within_distance(A_copy, B_copy, C_copy, D_copy, d);
}

am_interaction analyze_am(const vec& actin_left, const vec& actin_right,
    const vec& myosin_left, const vec& myosin_right,
    double& d, const std::vector<double>& box, const utils::PBCMask& pbc_mask)
{
    // Step 1: Copy segments
    vec actin_left_copy, actin_right_copy, myosin_left_copy, myosin_right_copy;
    copy_segments(actin_left, actin_right, myosin_left, myosin_right,
    actin_left_copy, actin_right_copy, myosin_left_copy, myosin_right_copy);

    // Step 2: Apply PBC on copies
    apply_pbc(actin_left_copy, actin_right_copy, myosin_left_copy, myosin_right_copy, box, pbc_mask);

    // Step 3: Analysis
    vec actin_mid = 0.5 * (actin_left_copy + actin_right_copy);
    vec myosin_mid = 0.5 * (myosin_left_copy + myosin_right_copy);

    std::tuple<double, vec, vec> result = subsegment_within_distance(
        actin_left_copy, actin_right_copy,
        myosin_left_copy, myosin_right_copy, d);

    am_interaction interaction;
    interaction.myosin_binding_ratio = std::get<0>(result);
    interaction.myosin_binding_start = std::get<1>(result);
    interaction.myosin_binding_end = std::get<2>(result);

    // Crosslinkable
    interaction.crosslinkable_start = actin_left_copy;
    interaction.crosslinkable_end = interaction.myosin_binding_start;
    interaction.crosslinkable_ratio = (interaction.crosslinkable_end - actin_left_copy).norm() /
                actin_right_copy.distance(actin_left_copy, box);

    if (interaction.crosslinkable_ratio < EPS) {
        interaction.partial_binding_ratio = 0.0;
        return interaction;
    }

    // Partial binding
    double dot = (actin_mid - actin_left_copy).dot(myosin_mid - myosin_left_copy);
    vec partial_start, partial_end;
    if (dot > 0) {
        partial_start = myosin_left_copy;
        partial_end = myosin_left_copy + (myosin_right_copy - myosin_left_copy) / 3.0;
    } else {
        partial_start = myosin_right_copy - (myosin_right_copy - myosin_left_copy) / 3.0;
        partial_end = myosin_right_copy;
    }

    std::tuple<double, vec, vec> result2 = subsegment_within_distance(
    actin_left_copy, actin_right_copy,
    partial_start, partial_end, d);

    interaction.partial_binding_ratio = std::get<0>(result2);

    if (interaction.partial_binding_ratio < EPS) {
        interaction.crosslinkable_ratio = 0.0;
    }

    return interaction;
}



} // namespace geometry
