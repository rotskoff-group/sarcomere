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

// Modified function to compute the shortest distance between a point and a segment (with PBC)
double point_segment_distance(vec x, vec a, vec b, std::vector<double> box)
{
    vec ab = b - a;
    double ab_norm = ab.norm();

    vec ab_normalized = {ab.x / ab_norm, ab.y / ab_norm};

    vec ap = x - a;
    ap.pbc_wrap(box);

    // Projection of ap onto ab
    double ap_ab = ap.x * ab_normalized.x + ap.y * ab_normalized.y;

    // Clamp projection length to segment bounds
    ap_ab = std::fmax(0, std::fmin(ap_ab, ab_norm));

    // Closest point on the segment
    vec closest_point = a + vec{ap_ab * ab_normalized.x, ap_ab * ab_normalized.y};

    // Vector from x to closest point
    vec norm_vec = x - closest_point;
    norm_vec.pbc_wrap(box);

    // Distance to the segment
    return norm_vec.norm();
}

std::pair<double,vec> point_segment_distance_w_normal(vec x, vec a, vec b, std::vector <double> box){
    // compute the shortest distance between a point x and a segment ab
    vec ab = b - a;
    double ab_norm = ab.norm();
    vec ab_normalized = {ab.x / ab_norm, ab.y / ab_norm};
    vec ap = x - a;
    ap.pbc_wrap(box);
    double ap_ab = ap.x * ab_normalized.x + ap.y * ab_normalized.y;
    ap_ab = std::fmax(0, std::fmin(ap_ab, ab_norm));
    vec norm_vec = ap - vec{ap_ab * ab_normalized.x, ap_ab * ab_normalized.y};
    double dist =  norm_vec.norm();
    return std::make_pair(dist, norm_vec);
}


// Compute orientation of three points with PBC
int orientation(const vec& p, const vec& q, const vec& r, const std::vector<double>& box)
{
    vec pq = q - p;
    vec qr = r - q;

    // Apply PBC wrapping
    pq.pbc_wrap(box);
    qr.pbc_wrap(box);

    // Calculate the orientation value
    double val = pq.y * qr.x - pq.x * qr.y;
    if (val == 0)
        return 0;  // Collinear
    return (val > 0) ? 1 : 2;  // Clockwise or Counterclockwise
}

// Compute orientation of three points without PBC
int orientation(const vec& p, const vec& q, const vec& r)
{
    vec pq = q - p;
    vec qr = r - q;

    // Calculate the orientation value
    double val = pq.y * qr.x - pq.x * qr.y;
    if (val == 0)
        return 0;  // Collinear
    return (val > 0) ? 1 : 2;  // Clockwise or Counterclockwise
}


// Compute the shortest distance between two segments with PBC
double segment_segment_distance(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box)
{
    // Segments are ab and cd
    double a_cd = point_segment_distance(a, c, d, box);
    double b_cd = point_segment_distance(b, c, d, box);
    double c_ab = point_segment_distance(c, a, b, box);
    double d_ab = point_segment_distance(d, a, b, box);

    double dist = std::fmin(std::fmin(a_cd, b_cd), std::fmin(c_ab, d_ab));

    // Check if segments intersect
    vec ab_center = {(a.x + b.x) / 2, (a.y + b.y) / 2};
    vec cd_center = {(c.x + d.x) / 2, (c.y + d.y) / 2};

    // Calculate displacement due to PBC
    vec displacement = {
        box[0] * std::round((ab_center.x - cd_center.x) / box[0]),
        box[1] * std::round((ab_center.y - cd_center.y) / box[1])
    };

    // Adjust a and b for displacement
    vec a0 = a - displacement;
    vec b0 = b - displacement;

    // Check orientations
    int o1 = orientation(a0, b0, c);
    int o2 = orientation(a0, b0, d);
    int o3 = orientation(c, d, a0);
    int o4 = orientation(c, d, b0);

    // If the orientations are such that the segments intersect, the distance is zero
    if ((o1 != o2) && (o3 != o4))
    {
        dist = 0;
    }

    return dist;
}


std::pair<double, std::map<std::string, vec>> segment_segment_distance_w_normal(const vec& a, const vec& b, const vec& c, const vec& d, const std::vector<double>& box) {
    // Compute point-to-segment distances with normal vectors
    std::pair<double, vec> a_cd = point_segment_distance_w_normal(a, c, d, box);
    std::pair<double, vec> b_cd = point_segment_distance_w_normal(b, c, d, box);
    std::pair<double, vec> c_ab = point_segment_distance_w_normal(c, a, b, box);
    std::pair<double, vec> d_ab = point_segment_distance_w_normal(d, a, b, box);

    // Get the distances
    double a_cd_dist = a_cd.first;
    double b_cd_dist = b_cd.first;
    double c_ab_dist = c_ab.first;
    double d_ab_dist = d_ab.first;

    // Find the minimum distance
    double dist = std::fmin(std::fmin(a_cd_dist, b_cd_dist), std::fmin(c_ab_dist, d_ab_dist));

    // Create a map to store the normal vector information
    std::map<std::string, vec> normal_dict;

    // Determine which point-to-segment distance is the smallest and store the corresponding normal vector
    if (dist == a_cd_dist) {
        normal_dict["start"] = a;
        normal_dict["vector"] = a_cd.second;
        normal_dict["end"] = a + a_cd.second;
    } else if (dist == b_cd_dist) {
        normal_dict["start"] = b;
        normal_dict["vector"] = b_cd.second;
        normal_dict["end"] = b + b_cd.second;
    } else if (dist == c_ab_dist) {
        normal_dict["end"] = c;
        normal_dict["vector"] = -c_ab.second;
        normal_dict["start"] = c - c_ab.second;
    } else {
        normal_dict["end"] = d;
        normal_dict["vector"] = -d_ab.second;
        normal_dict["start"] = d - d_ab.second;
    }

    // Check if the segments intersect
    vec ab_center = {(a.x + b.x) / 2, (a.y + b.y) / 2};
    vec cd_center = {(c.x + d.x) / 2, (c.y + d.y) / 2};

    // Calculate displacement due to PBC
    vec displacement = {
        box[0] * std::round((ab_center.x - cd_center.x) / box[0]),
        box[1] * std::round((ab_center.y - cd_center.y) / box[1])
    };

    // Adjust a and b for displacement
    vec a0 = a - displacement;
    vec b0 = b - displacement;

    // Check orientations to determine if the segments intersect
    int o1 = orientation(a0, b0, c, box);
    int o2 = orientation(a0, b0, d, box);
    int o3 = orientation(c, d, a0, box);
    int o4 = orientation(c, d, b0, box);

    // If the segments intersect, set the distance to 0 and the normal vector to zero
    if ((o1 != o2) && (o3 != o4)) {
        dist = 0;
        normal_dict["vector"] = {0, 0};
    }

    return std::make_pair(dist, normal_dict);

}
// Function to compute the dot product of two vectors
double dotProduct(double x1, double y1, double x2, double y2) {
    return x1 * x2 + y1 * y2;
}



void apply_pbc(vec& actin_left, vec& actin_right, 
                        vec& myosin_left, vec& myosin_right, std::vector<double> box){

    vec actin_mid = (actin_left + actin_right) / 2;
    vec myosin_mid = (myosin_left + myosin_right) / 2;
    vec displacement = actin_mid - myosin_mid;
    vec shift;
    shift.x =  -box[0]*round(displacement.x/box[0]);
    shift.y =  - box[1]*round(displacement.y/box[1]);
    actin_left = actin_left + shift;
    actin_right = actin_right + shift;
}




// Filter roots within the interval [t_start, t_end]
std::vector<double> findRootsInInterval(const std::vector<double>& roots, double t_start, double t_end) {
    std::vector<double> valid_roots;
    for (double t : roots) {
        if (t >= t_start - EPS && t <= t_end + EPS) {
            valid_roots.push_back(t);
        }
    }
    return valid_roots;
}




// Solve quadratic equation: a*t^2 + b*t + c = 0
std::vector<double> solveQuadratic(double a, double b, double c) {
    std::vector<double> roots;
    double discriminant = b * b - 4 * a * c;

    if (discriminant < -EPS) {
        // No real roots
        return roots;
    } else if (std::abs(discriminant) < EPS) {
        // One real root
        double t = -b / (2 * a);
        roots.push_back(t);
    } else {
        // Two real roots
        double sqrt_disc = std::sqrt(discriminant);
        double t1 = (-b + sqrt_disc) / (2 * a);
        double t2 = (-b - sqrt_disc) / (2 * a);
        roots.push_back(t1);
        roots.push_back(t2);
    }
    return roots;
}


// Main function to find intervals on AB where distance to CD is less than d
std::tuple<double,vec,vec> subsegment_within_distance(vec A, vec B, vec C, vec D, double d) {
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
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
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
            for (double t_root : roots) {
                if (t_root >= t_start - EPS && t_root <= t_end + EPS) {
                    valid_roots.push_back(t_root);
                }
            }

            // Collect interval points
            std::vector<double> interval_points = {t_start};
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
    // // Merge overlapping intervals
    // if (!intervals.empty()) {
    //     // Sort intervals by t_start
    //     std::sort(intervals.begin(), intervals.end(), [&](const std::pair<vec, vec>& a, const std::pair<vec, vec>& b) {
    //         double t_a = (a.first - A).dot(dirAB) / lenAB_sq;
    //         double t_b = (b.first - A).dot(dirAB) / lenAB_sq;
    //         return t_a < t_b;
    //     });

    //     std::vector<std::pair<vec, vec>> merged_intervals;
    //     auto current = intervals[0];

    //     for (size_t i = 1; i < intervals.size(); ++i) {
    //         double t_current_end = (current.second - A).dot(dirAB) / lenAB_sq;
    //         double t_next_start = (intervals[i].first - A).dot(dirAB) / lenAB_sq;

    //         if (t_next_start <= t_current_end + EPS) {
    //             // Overlapping intervals, merge them
    //             current.second = intervals[i].second;
    //         } else {
    //             merged_intervals.push_back(current);
    //             current = intervals[i];
    //         }
    //     }
    //     merged_intervals.push_back(current);
    //     return merged_intervals;
    // }

    // return intervals;
}



std::tuple<double,vec,vec>subsegment_within_distance(vec A, vec B, vec C, vec D, double d, std::vector<double> box) {
    apply_pbc(A, B, C, D, box);
    return subsegment_within_distance(A, B, C, D, d);
}


// Function to calculate the points on segment 1 that are distance d away from segment 2
am_interaction analyze_am(vec actin_left, vec actin_right, 
                        vec myosin_left, vec myosin_right, double d, std::vector<double> box) {
    apply_pbc(actin_left, actin_right, myosin_left, myosin_right, box);
    vec actin_mid = 0.5 * (actin_left + actin_right);
    vec myosin_mid = 0.5 * (myosin_left + myosin_right);
    std::tuple<double, vec, vec> result = subsegment_within_distance(actin_left, actin_right, myosin_left, myosin_right, d);
    am_interaction interaction;
    interaction.myosin_binding_ratio = std::get<0>(result);
    interaction.myosin_binding_start = std::get<1>(result);
    interaction.myosin_binding_end = std::get<2>(result);

    //crosslinkable
    interaction.crosslinkable_start = actin_left;
    interaction.crosslinkable_end = interaction.myosin_binding_start;
    interaction.crosslinkable_ratio = (interaction.crosslinkable_end-actin_left).norm()/actin_right.distance(actin_left,box);
    if (interaction.crosslinkable_ratio<EPS){
        interaction.partial_binding_ratio = 0.0;
        return interaction;
    }
    double dot = (actin_mid-actin_left).dot(myosin_mid-myosin_left);
    vec partial_start, partial_end;
    if (dot > 0){
        partial_start = myosin_left;
        partial_end = myosin_left + (myosin_right-myosin_left)/3.0;
    }
    else{
        partial_start = myosin_right - (myosin_right-myosin_left)/3.0;
        partial_end = myosin_right;
    }
    std::tuple<double, vec, vec> result2 = subsegment_within_distance(actin_left, actin_right, partial_start, partial_end, d);
    interaction.partial_binding_ratio = std::get<0>(result2);
    if (interaction.partial_binding_ratio < EPS){
        interaction.crosslinkable_ratio = 0.0;
    }
    return interaction;   
}

} // namespace geometry








