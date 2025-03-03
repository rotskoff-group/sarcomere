#include "interaction.h"
#include <cstdio>
#include <cstdlib>  // For exit()

//---------------------------------------------------------------------
// Non-Template Function Definitions
//---------------------------------------------------------------------

real aa_energy(const ArrayXreal& center1, const double& length1, const real& theta1,
               const ArrayXreal& center2, const double& length2, const real& theta2,
               const std::vector<double>& box, const double k_aa, const double kappa_aa)
{
    // Calculate the endpoints of the first segment.
    real x1_start = center1[0] - (length1 / 2) * cos(theta1);
    real y1_start = center1[1] - (length1 / 2) * sin(theta1);
    real x1_end   = center1[0] + (length1 / 2) * cos(theta1);
    real y1_end   = center1[1] + (length1 / 2) * sin(theta1);

    // Calculate the endpoints of the second segment.
    real x2_start = center2[0] - (length2 / 2) * cos(theta2);
    real y2_start = center2[1] - (length2 / 2) * sin(theta2);
    real x2_end   = center2[0] + (length2 / 2) * cos(theta2);
    real y2_end   = center2[1] + (length2 / 2) * sin(theta2);

    real a[2] = {x1_start, y1_start};
    real b[2] = {x1_end,   y1_end};
    real c[2] = {x2_start, y2_start};
    real d[2] = {x2_end,   y2_end};

    // Compute distance between segments.
    real dist = segment_segment_distance(a, b, c, d, box);
    dist = dist - 0.03;
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = std::min(abs(angle), M_PI - abs(angle));
    return 0.5 * (k_aa * dist * dist + kappa_aa * angle * angle);
}

real am_energy1(const ArrayXreal& center1, const double& length1, const real& theta1,
                const ArrayXreal& center2, const double& length2, const real& theta2,
                const std::vector<double>& box, const double k_am, const double kappa_am,
                const double myosin_radius)
{
    // Calculate the endpoints of the first segment.
    real x1_start = center1[0] - (length1 / 2) * cos(theta1);
    real y1_start = center1[1] - (length1 / 2) * sin(theta1);
    real x1_end   = center1[0] + (length1 / 2) * cos(theta1);
    real y1_end   = center1[1] + (length1 / 2) * sin(theta1);

    // Calculate the endpoints of the second segment.
    real x2_start = center2[0] - (length2 / 2) * cos(theta2);
    real y2_start = center2[1] - (length2 / 2) * sin(theta2);
    real x2_end   = center2[0] + (length2 / 2) * cos(theta2);
    real y2_end   = center2[1] + (length2 / 2) * sin(theta2);

    real a[2] = {x1_start, y1_start};
    real b[2] = {x1_end,   y1_end};
    real c[2] = {x2_start, y2_start};
    real d[2] = {x2_end,   y2_end};

    real dist = segment_segment_distance(a, b, c, d, box);
    real strength = abs(cos(theta1 - theta2));
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = abs(angle);
    angle = (angle < M_PI - angle) ? angle : (M_PI - angle);
    if (dist.val() > 0.8 * myosin_radius) {
        if (dist.val() > myosin_radius) {
            printf("dist: %f\n", dist.val());
            printf("endpoints after: %f %f %f %f %f %f %f %f\n",
                   val(a[0]), val(a[1]), val(b[0]), val(b[1]),
                   val(c[0]), val(c[1]), val(d[0]), val(d[1]));
            printf("something is wrong\n");
        }
        return 0.5 * k_am * val(strength) * dist * dist + 0.5 * kappa_am * angle * angle;
    }
    else {
        real energy = 0.5 * (k_am / 10) * val(strength) * dist * dist + 0.5 * kappa_am * angle * angle;
        return energy;
    }
}

real am_energy(const real& theta1, const real& theta2, const double kappa_am)
{
    real angle = theta1 - theta2;
    angle_wrap(angle);
    angle = abs(angle);
    angle = (angle < M_PI - angle) ? angle : (M_PI - angle);
    return 0.5 * kappa_am * angle * angle;
}

std::vector<double> compute_aa_force_and_energy(Filament& actin,
                                                int& actin1_index, int& actin2_index,
                                                const std::vector<double>& box,
                                                const double k_aa, const double kappa_aa)
{
    ArrayXreal center1(2);
    center1 << actin.center[actin1_index].x, actin.center[actin1_index].y;
    real theta1 = actin.theta[actin1_index];
    ArrayXreal center2(2);
    center2 << actin.center[actin2_index].x, actin.center[actin2_index].y;
    real theta2 = actin.theta[actin2_index];
    real u;
    std::vector<double> forces;
    VectorXd forces_2 = -gradient(aa_energy, wrt(center1, theta1, theta2),
                                    at(center1, actin.length, theta1,
                                       center2, actin.length, theta2, box, k_aa, kappa_aa), u);
    forces.resize(forces_2.size());
    VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    return forces;
}


std::vector<double> compute_am_force_and_energy(Filament& actin, Myosin& myosin,
                                                int& actin_index, int& myosin_index,
                                                const std::vector<double>& box,
                                                const double k_am, const double kappa_am,
                                                const double myosin_radius)
{
    ArrayXreal center1(2);
    center1 << actin.center[actin_index].x, actin.center[actin_index].y;
    real theta1 = actin.theta[actin_index];
    ArrayXreal center2(2);
    center2 << myosin.center[myosin_index].x, myosin.center[myosin_index].y;
    real theta2 = myosin.theta[myosin_index];
    real u;
    std::vector<double> forces;
    VectorXd forces_2;
    if (k_am > 1e-6) {
        forces_2 = -gradient(am_energy1, wrt(center1, theta1, theta2),
                             at(center1, actin.length, theta1,
                                center2, myosin.length, theta2, box, k_am, kappa_am, myosin_radius), u);
        forces.resize(forces_2.size());
        // (Energy value computed here is not used in the force vector.)
        VectorXd::Map(&forces[0], forces_2.size()) = forces_2;
    }
    else {
        forces_2 = -gradient(am_energy, wrt(theta1, theta2),
                             at(theta1, theta2, kappa_am), u);
        // Prepend two zeros (if needed) for the force vector.
        forces.resize(forces_2.size() + 2);
        forces[0] = 0;
        forces[1] = 0;
        VectorXd::Map(&forces[2], forces_2.size()) = forces_2;
    }
    return forces;
}