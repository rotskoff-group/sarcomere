#include "interaction.h"
#include <cmath>
#include <vector>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <Eigen/Dense>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>

using namespace autodiff;
using Eigen::VectorXd;
using std::vector;

real smooth_round(real x) {
    return round(val(x));
}

void angle_wrap(real& theta) {
    theta = theta - 2.0 * M_PI * smooth_round(theta / (2.0 * M_PI));
}


vector<real> pbc_wrap(const vector<real>& x, const vector<double>& box) {
    vector<real> wrapped(2);
    wrapped[0] = x[0] - box[0] * smooth_round(x[0] / box[0]);
    wrapped[1] = x[1] - box[1] * smooth_round(x[1] / box[1]);
    return wrapped;
}


void apply_pbc(vector<real>& a, vector<real>& b, vector<real>& c, vector<real>& d,
               const vector<double>& box, const utils::PBCMask& pbc_mask) {
    vector<real> ab_center = {(a[0] + b[0]) / 2.0, (a[1] + b[1]) / 2.0};
    vector<real> cd_center = {(c[0] + d[0]) / 2.0, (c[1] + d[1]) / 2.0};
    real shift;
    if (pbc_mask.x) {
        real dx = ab_center[0] - cd_center[0];
        shift = box[0] * smooth_round(dx / box[0]);
        a[0] -= shift; b[0] -= shift;
    }
    if (pbc_mask.y) {
        real dy = ab_center[1] - cd_center[1];
        shift = box[1] * smooth_round(dy / box[1]);
        a[1] -= shift; b[1] -= shift;
    }
}

real point_segment_distance(const vector<real>& x, const vector<real>& a, const vector<real>& b,
                             const vector<double>& box, vector<real>& norm_vec) {
    vector<real> ab = {b[0] - a[0], b[1] - a[1]};
    real ab_norm = sqrt(ab[0]*ab[0] + ab[1]*ab[1]);
    vector<real> ab_normalized = {ab[0] / ab_norm, ab[1] / ab_norm};

    vector<real> ap = {x[0] - a[0], x[1] - a[1]};
    ap = pbc_wrap(ap, box);

    real ap_ab = ap[0] * ab_normalized[0] + ap[1] * ab_normalized[1];
    ap_ab = max(real(0), min(ap_ab, ab_norm));

    norm_vec = {ap[0] - ap_ab * ab_normalized[0],
                ap[1] - ap_ab * ab_normalized[1]};
    return sqrt(norm_vec[0] * norm_vec[0] + norm_vec[1] * norm_vec[1]);
}

real point_segment_distance(const vector<real>& x, const vector<real>& a, const vector<real>& b,
                             vector<real>& norm_vec) {
    vector<real> ab = {b[0] - a[0], b[1] - a[1]};
    real ab_norm = sqrt(ab[0]*ab[0] + ab[1]*ab[1]);
    vector<real> ab_normalized = {ab[0] / ab_norm, ab[1] / ab_norm};

    vector<real> ap = {x[0] - a[0], x[1] - a[1]};

    real ap_ab = ap[0] * ab_normalized[0] + ap[1] * ab_normalized[1];
    ap_ab = max(real(0), min(ap_ab, ab_norm));

    norm_vec = {ap[0] - ap_ab * ab_normalized[0],
                ap[1] - ap_ab * ab_normalized[1]};
    return sqrt(norm_vec[0] * norm_vec[0] + norm_vec[1] * norm_vec[1]);
}

int orientation(const vector<real>& p, const vector<real>& q, const vector<real>& r) {
    real pq_x = q[0] - p[0];
    real pq_y = q[1] - p[1];
    real qr_x = r[0] - q[0];
    real qr_y = r[1] - q[1];
    real val = pq_y * qr_x - pq_x * qr_y;
    if (val == real(0)) return 0;
    return (val > 0) ? 1 : 2;
}

real segment_segment_distance(const vector<real>& a_in, const vector<real>& b_in,
                               const vector<real>& c_in, const vector<real>& d_in,
                               const vector<double>& box, const utils::PBCMask& pbc_mask) {
    vector<real> a = a_in, b = b_in, c = c_in, d = d_in;
    
    apply_pbc(a, b, c, d, box, pbc_mask);
    vector<real> norm_vec;
    real a_cd = point_segment_distance(a, c, d, norm_vec);
    real b_cd = point_segment_distance(b, c, d, norm_vec);
    real c_ab = point_segment_distance(c, a, b, norm_vec);
    real d_ab = point_segment_distance(d, a, b, norm_vec);
    real dist = min(min(a_cd, b_cd), min(c_ab, d_ab));

    int o1 = orientation(a, b, c);
    int o2 = orientation(a, b, d);
    int o3 = orientation(c, d, a);
    int o4 = orientation(c, d, b);

    if ((o1 != o2) && (o3 != o4)) {
        dist = real(0);
    }

    return dist;
}

real aa_energy(const ArrayXreal& center1, const double& length1, const real& theta1,
               const ArrayXreal& center2, const double& length2, const real& theta2,
               const vector<double>& box, const utils::PBCMask& pbc_mask,
               const double k_aa, const double kappa_aa)
{
    vector<real> a = {center1[0] - (length1/2) * cos(theta1),
                      center1[1] - (length1/2) * sin(theta1)};
    vector<real> b = {center1[0] + (length1/2) * cos(theta1),
                      center1[1] + (length1/2) * sin(theta1)};
    vector<real> c = {center2[0] - (length2/2) * cos(theta2),
                      center2[1] - (length2/2) * sin(theta2)};
    vector<real> d = {center2[0] + (length2/2) * cos(theta2),
                      center2[1] + (length2/2) * sin(theta2)};

    real dist = segment_segment_distance(a, b, c, d, box, pbc_mask);
    dist = dist - 0.03;
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = min(abs(angle), M_PI - abs(angle));
    return 0.5 * (k_aa * dist * dist + kappa_aa * angle * angle);
}

real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1,
                const ArrayXreal& center2, const double& length2, const real& theta2,
                const vector<double>& box, const utils::PBCMask& pbc_mask,
                const double k_am, const double kappa_am,
                const double myosin_radius)
{
    vector<real> a = {center1[0] - (length1/2) * cos(theta1),
                      center1[1] - (length1/2) * sin(theta1)};
    vector<real> b = {center1[0] + (length1/2) * cos(theta1),
                      center1[1] + (length1/2) * sin(theta1)};
    vector<real> c = {center2[0] - (length2/2) * cos(theta2),
                      center2[1] - (length2/2) * sin(theta2)};
    vector<real> d = {center2[0] + (length2/2) * cos(theta2),
                      center2[1] + (length2/2) * sin(theta2)};

    real dist = segment_segment_distance(a, b, c, d, box, pbc_mask);
    
    real strength = abs(cos(theta1 - theta2));
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = min(abs(angle), M_PI - abs(angle));
    if (dist.val() > 0.85 * myosin_radius) {
        real offset = dist - 0.85 * myosin_radius;
        return 0.5 * k_am * val(strength) * offset * offset + 0.5 * kappa_am * angle * angle;
    } else {
        return 0.5 * kappa_am * angle * angle;
    }
}

real am_energy(const real& theta1, const real& theta2, const double kappa_am)
{
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = min(abs(angle), M_PI - abs(angle));
    return 0.5 * kappa_am * angle * angle;
}

vector<double> compute_aa_force_and_energy(Filament& actin, int& i, int& j,
                                           const vector<double>& box, const utils::PBCMask& pbc_mask,
                                           const double k_aa, const double kappa_aa)
{
    ArrayXreal center1(2);
    center1 << actin.center[i].x, actin.center[i].y;
    real theta1 = actin.theta[i];
    ArrayXreal center2(2);
    center2 << actin.center[j].x, actin.center[j].y;
    real theta2 = actin.theta[j];
    real energy;
    VectorXd grad = -gradient(aa_energy, wrt(center1, theta1, theta2),
                              at(center1, actin.length, theta1, center2, actin.length, theta2, box, pbc_mask, k_aa, kappa_aa), energy);
    vector<double> forces(grad.size());
    Eigen::Map<VectorXd>(&forces[0], grad.size()) = grad;
    return forces;
}

vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin, int& i, int& j,
                                           const vector<double>& box, const utils::PBCMask& pbc_mask,
                                           const double k_am, const double kappa_am, const double myosin_radius)
{
    ArrayXreal center1(2);
    center1 << actin.center[i].x, actin.center[i].y;
    real theta1 = actin.theta[i];
    ArrayXreal center2(2);
    center2 << myosin.center[j].x, myosin.center[j].y;
    real theta2 = myosin.theta[j];
    real energy;
    VectorXd grad;

    if (k_am > 1e-6) {
        grad = -gradient(am_energy1, wrt(center1, theta1, theta2),
                         at(center1, actin.length, theta1,
                            center2, myosin.length, theta2, box, pbc_mask, k_am, kappa_am, myosin_radius), energy);
        vector<double> forces(grad.size());
        Eigen::Map<VectorXd>(&forces[0], grad.size()) = grad;
        //check if the force is too large
        // if (forces[0] > 3 || forces[1] > 3) {
        //     real energy = am_energy1(center1, actin.length, theta1,
        //         center2, myosin.length, theta2, box, pbc_mask, k_am, kappa_am, myosin_radius);
        //     std::cout << "am force too large: " << forces[0] << " " << forces[1] << " energy: " << energy << std::endl;
        // }
        return forces;
    } else {
        grad = -gradient(am_energy, wrt(theta1, theta2), at(theta1, theta2, kappa_am), energy);
        vector<double> forces(2 + grad.size(), 0.0);
        Eigen::Map<VectorXd>(&forces[2], grad.size()) = grad;
        return forces;
    }
}
