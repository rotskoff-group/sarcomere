#include "components.h"
#include <cstdio>
#include <stdexcept>
#include <cmath>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

//===================
// Filament Methods
//===================

// Default constructor.
Filament::Filament() : n(0), length(0) {
    // Optionally initialize other members if needed.
}

// Parameterized constructor.
Filament::Filament(int n0, double length0, std::vector<double> box0, utils::PBCMask pbc_mask0, gsl_rng* rng)
    : n(n0), length(length0), box(box0), pbc_mask(pbc_mask0)
{
    center.resize(n);
    theta.resize(n);
    left_end.resize(n);
    right_end.resize(n);
    force.resize(n);
    angular_force.resize(n);
    velocity.resize(n);
    f_load.resize(n);
    cb_strength.resize(n);
    
    // Randomly initialize the center positions and theta values.
    for (int i = 0; i < n; i++) {
        theta[i] = gsl_ran_flat(rng, 0, 2 * M_PI);
        if (!pbc_mask.x) {
            double dx = 0.5 * length * std::abs(std::cos(theta[i]));
            double x_min = -0.5 * box[0] + dx;
            double x_max =  0.5 * box[0] - dx;
            center[i].x = gsl_ran_flat(rng, x_min, x_max);
        } else {
            center[i].x = gsl_ran_flat(rng, -0.5 * box[0], 0.5 * box[0]);
        }
        if (!pbc_mask.y) {
            double dy = 0.5 * length * std::abs(std::sin(theta[i]));
            double y_min = -0.5 * box[1] + dy;
            double y_max =  0.5 * box[1] - dy;
            center[i].y = gsl_ran_flat(rng, y_min, y_max);
        } else {
            center[i].y = gsl_ran_flat(rng, -0.5 * box[1], 0.5 * box[1]);
        }
    }

    update_endpoints();
}

// Destructor.
Filament::~Filament() {
    printf("Filament destructor called\n");
}

// Copy constructor.
Filament::Filament(const Filament& other) {
    n = other.n;
    box = other.box;
    length = other.length;
    cb_strength = other.cb_strength;
    center = other.center;
    theta = other.theta;
    left_end = other.left_end;
    right_end = other.right_end;
    force = other.force;
    angular_force = other.angular_force;
    velocity = other.velocity;
    f_load = other.f_load;
    custom_features = other.custom_features;
    pbc_mask = other.pbc_mask;
}

void Filament::displace(int& i, double& dx, double& dy, double& dtheta) {
    if (pbc_mask.x && pbc_mask.y) {
        center[i].x += dx;
        center[i].y += dy;
        center[i].pbc_wrap(box);
        theta[i] += dtheta;
        utils::angle_wrap(theta[i]);
        update_endpoints(i);
    } else{
        // Proposed new center and angle
        vec proposed_center = center[i];
        proposed_center.x += dx;
        proposed_center.y += dy;

        double proposed_theta = theta[i] + dtheta;
        utils::angle_wrap(proposed_theta);

        // Compute proposed endpoints
        double dx_seg = 0.5 * length * std::cos(proposed_theta);
        double dy_seg = 0.5 * length * std::sin(proposed_theta);
        double y1 = proposed_center.y - dy_seg;
        double y2 = proposed_center.y + dy_seg;
        if (!pbc_mask.x) {
            double x1 = proposed_center.x - dx_seg;
            double x2 = proposed_center.x + dx_seg;
            // Check if either endpoint is outside the x-boundaries
            if (x1 < -0.5 * box[0] || x1 > 0.5 * box[0] ||
                x2 < -0.5 * box[0] || x2 > 0.5 * box[0]) {
                return; // Reject move
            }
        }
        if (!pbc_mask.y) {
            // Check if either endpoint is outside the y-boundaries
            if (y1 < -0.5 * box[1] || y1 > 0.5 * box[1] ||
                y2 < -0.5 * box[1] || y2 > 0.5 * box[1]) {
                return; // Reject move
            }
        }
        // Accept move
        center[i] = proposed_center;
        center[i].pbc_wrap(box, pbc_mask);
        theta[i] = proposed_theta;
        update_endpoints(i);
    }
}


// Update endpoints for the i-th filament.
void Filament::update_endpoints(int& i) {
    std::vector<double> segments(2);
    segments[0] = length * cos(theta[i]);
    segments[1] = length * sin(theta[i]);
    left_end[i].x = center[i].x - 0.5 * segments[0];
    left_end[i].y = center[i].y - 0.5 * segments[1];
    right_end[i].x = center[i].x + 0.5 * segments[0];
    right_end[i].y = center[i].y + 0.5 * segments[1];
}

// Update endpoints for all filaments.
void Filament::update_endpoints() {
    for (int i = 0; i < n; i++) {
        update_endpoints(i);
    }
}

// Update theta for all filaments.
void Filament::update_theta(std::vector<double> new_theta) {
    for (int i = 0; i < n; i++) {
        theta[i] = new_theta[i];
    }
    update_endpoints();
}

// Update center positions for all filaments.
void Filament::update_center(std::vector<vec> new_center) {
    for (int i = 0; i < n; i++) {
        center[i] = new_center[i];
    }
    update_endpoints();
}

// Register a new 1D feature.
void Filament::register_feature(const std::string& name) {
    if (custom_features.find(name) == custom_features.end()) {
        custom_features[name] = std::vector<double>(n, 0.0); // Initialize with zeros.
        std::cout << "1D feature " << name << " registered successfully.\n";
    } else {
        std::cout << "1D feature " << name << " already exists.\n";
    }
}

// Overload operator[] to access custom features.
std::vector<double>& Filament::operator[](const std::string& name) {
    if (custom_features.find(name) == custom_features.end()) {
        throw std::runtime_error("Error: 1D feature " + name + " not found.");
    }
    return custom_features[name];
}

//===================
// Myosin Methods
//===================

// Default constructor.
Myosin::Myosin() : Filament() {
    // The base default constructor is automatically called.
}

// Parameterized constructor.
Myosin::Myosin(int n0, double length0, double radius0, std::vector<double> box0, utils::PBCMask pbc_mask0, gsl_rng* rng)
    : Filament(n0, length0, box0, pbc_mask0, rng), radius(radius0)
{
    // Additional initialization if needed.
}

// Copy constructor.
Myosin::Myosin(const Myosin& other) : Filament(other) {
    radius = other.radius;
}
