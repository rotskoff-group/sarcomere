#ifndef COMPONENTS_H
#define COMPONENTS_H

#ifndef UTILS_H
#include "utils.h"
#endif

#include <cmath>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <string>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

// Define an alias for vec using the type from utils.
using vec = utils::vec;

class Filament {
public:
    int n;
    std::vector<double> box;
    double length;
    std::vector<vec> center;
    std::vector<double> theta;
    std::vector<vec> left_end;
    std::vector<vec> right_end;
    std::vector<vec> force;
    std::vector<double> angular_force;
    std::vector<vec> velocity;
    std::vector<double> f_load;
    std::vector<double> cb_strength;

    // Dictionary to store 1D features.
    std::unordered_map<std::string, std::vector<double>> custom_features;

    // Constructors and destructor.
    Filament();
    Filament(int n0, double length0, std::vector<double> box0, gsl_rng* rng);
    virtual ~Filament();
    Filament(const Filament& other);

    // Member functions.
    void displace(int& i, double& dx, double& dy);
    void displace(int& i, double& dx, double& dy, double& dtheta);
    void update_endpoints(int& i);
    void update_endpoints();
    void update_theta(std::vector<double> new_theta);
    void update_center(std::vector<vec> new_center);

    // Register a new 1D feature of length n.
    void register_feature(const std::string& name);

    // Operator overload for easy access to 1D features.
    std::vector<double>& operator[](const std::string& name);
};

class Myosin : public Filament {
public:
    double radius;

    // Constructors.
    Myosin();
    Myosin(int n0, double length0, double radius0, std::vector<double> box0, gsl_rng* rng);
    Myosin(const Myosin& other);
};

#endif // COMPONENTS_H
