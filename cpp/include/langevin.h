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
    Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin, double& D0_myosin,
             int& update_myosin_every0, int& update_dt_every0, int& save_every0, bool& resume);

    // Destructor.
    ~Langevin();

    // Run the Langevin simulation for a given number of steps.
    void run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin);
    void volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin);
    // Take a single simulation step.
    void sample_step(double& dt, double& D_actin, double& D_myosin, gsl_rng* rng, int& fix_myosin);

    // Data members.
    double dt, beta, acc_rate, D_actin, D_myosin;
    int update_myosin_every, update_dt_every, save_every;
    Sarcomere& model;
};

#endif // LANGEVIN_H
