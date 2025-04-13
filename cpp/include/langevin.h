#ifndef LANGEVIN_H
#define LANGEVIN_H

#include "sarcomere.h"
#include "utils.h"
#include "h5_utils.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <iostream>
#include <omp.h>
#include <mutex>

// Declare the global mutex for saving (defined in the source file)
extern std::mutex save_mutex;

class Langevin {
public:
    // Constructor.
    // Note: Parameters are passed by reference as in the original code.
    Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
        double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot,
        int& save_every0, bool& resume);
    // Destructor.
    ~Langevin();

    // Run the Langevin simulation for a given number of steps.
    void run_langevin(int& nsteps, gsl_rng* rng, int& fix_myosin);
    void volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin);
    // Take a single simulation step.
    void sample_step(double& dt, gsl_rng* rng, int& fix_myosin);

    // Data members.
    double dt, beta, acc_rate, D_actin_trans, D_actin_rot, D_myosin_trans, D_myosin_rot;
    int save_every, start_step;
    Sarcomere& model;
};

#endif // LANGEVIN_H
