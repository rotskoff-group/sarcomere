#ifndef INTERACTION_H
#define INTERACTION_H

#include <vector>
#include <autodiff/forward/real.hpp>
#include <autodiff/forward/real/eigen.hpp>
#include "components.h"
#include "utils.h"  // Defines utils::PBCMask

using autodiff::real;
using autodiff::ArrayXreal;
using std::vector;

//
// Function declarations
//

// Math utilities
real smooth_round(real x);
void angle_wrap(real& theta);
vector<real> pbc_wrap(const vector<real>& x, const vector<double>& box);

// Periodic boundary handling
void apply_pbc(vector<real>& a, vector<real>& b, vector<real>& c, vector<real>& d,
               const vector<double>& box, const utils::PBCMask& pbc_mask);

// Geometric distances
real point_segment_distance(const vector<real>& x, const vector<real>& a, const vector<real>& b,
                             const vector<double>& box, vector<real>& norm_vec);
real point_segment_distance(const vector<real>& x, const vector<real>& a, const vector<real>& b);

// Orientation
int orientation(const vector<real>& p, const vector<real>& q, const vector<real>& r);

// Segment–segment distance
real segment_segment_distance(const vector<real>& a, const vector<real>& b,
                               const vector<real>& c, const vector<real>& d,
                               const vector<double>& box, const utils::PBCMask& pbc_mask);

// Energies
real aa_energy(const ArrayXreal& center1, const double& length1, const real& theta1,
               const ArrayXreal& center2, const double& length2, const real& theta2,
               const vector<double>& box, const utils::PBCMask& pbc_mask,
               const double k_aa, const double kappa_aa);

real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1,
                const ArrayXreal& center2, const double& length2, const real& theta2,
                const vector<double>& box, const utils::PBCMask& pbc_mask,
                const double k_am, const double kappa_am, const double myosin_radius);

real am_energy(const real& theta1, const real& theta2, const double kappa_am);

// Force computations
vector<double> compute_aa_force_and_energy(Filament& actin, int& i, int& j,
                                           const vector<double>& box, const utils::PBCMask& pbc_mask,
                                           const double k_aa, const double kappa_aa);

vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin, int& i, int& j,
                                           const vector<double>& box, const utils::PBCMask& pbc_mask,
                                           const double k_am, const double kappa_am,
                                           const double myosin_radius);

#endif  // INTERACTION_H
