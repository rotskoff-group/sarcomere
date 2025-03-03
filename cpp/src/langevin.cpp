#include "langevin.h"
#include <cstdio>
#include <cmath>
#include <vector>

// Define the global mutex.
std::mutex save_mutex;

//---------------------------------------------------------------------
// Constructor: Initializes the Langevin object and, if resume is true,
// loads the previous state and prints details. Otherwise, creates a new file.
//---------------------------------------------------------------------
Langevin::Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin, double& D0_myosin,
                   int& update_myosin_every0, int& update_dt_every0, int& save_every0, bool& resume)
    : model(model0), beta(beta0), dt(dt0), D_actin(D0_actin), D_myosin(D0_myosin),
      update_myosin_every(update_myosin_every0), update_dt_every(update_dt_every0),
      save_every(save_every0), acc_rate(0)
{
    if (resume) {
        model.load_state();
        printf("Resuming from file %s\n", model.filename.c_str());
        // Print myosin details.
        for (int i = 0; i < model.myosin.n; i++) {
            printf("Myosin %d: %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y);
            printf("Myosin endpoints: (%f %f), (%f %f)\n",
                   model.myosin.left_end[i].x, model.myosin.left_end[i].y,
                   model.myosin.right_end[i].x, model.myosin.right_end[i].y);
        }
        // Print actin details.
        for (int i = 0; i < model.actin.n; i++) {
            printf("Actin %d: %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y);
            printf("Actin endpoints: (%f %f), (%f %f)\n",
                   model.actin.left_end[i].x, model.actin.left_end[i].y,
                   model.actin.right_end[i].x, model.actin.right_end[i].y);
        }
    } else {
        model.new_file();
    }
}

//---------------------------------------------------------------------
// Destructor
//---------------------------------------------------------------------
Langevin::~Langevin() {
    // No dynamic resources need explicit deallocation.
}

//---------------------------------------------------------------------
// run_langevin: Runs the simulation for nsteps, periodically saving the
// state and taking sample steps.
//---------------------------------------------------------------------
void Langevin::run_langevin(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    for (int i = 0; i < nsteps; i++) {
        if (i % save_every == 0) {
            std::cout << "Step " << i << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            model.save_state();
            start = omp_get_wtime();
        }
        model.update_system();
        sample_step(dt, D_actin, D_myosin, rng, fix_myosin);
        if (i % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", i, end - start);
        }
    }
}

void Langevin::volume_exclusion(int nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    for (int i = 0; i < nsteps; i++) {
        if (i % save_every == 0) {
            std::cout << "Step " << i << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            start = omp_get_wtime();
        }
        model.update_system_sterics_only();
        sample_step(dt, D_actin, D_myosin, rng, fix_myosin);
        if (i % save_every == 0) {
            end = omp_get_wtime();
            printf("Step %d took %f seconds\n", i, end - start);
        }
    }
}

//---------------------------------------------------------------------
// sample_step: Performs a single Langevin dynamics step by updating the 
// system, generating noise, and displacing myosin and actin particles.
//---------------------------------------------------------------------
void Langevin::sample_step(double& dt, double& D_actin, double& D_myosin, gsl_rng* rng, int& fix_myosin) {
    // Generate noise for both myosin and actin particles.
    int n_randns = (model.myosin.n + model.actin.n) * 3;
    std::vector<double> noise(n_randns);
    for (int i = 0; i < n_randns; i++) {
        noise[i] = gsl_rng_uniform(rng) * 2 - 1;
    }

    int n_acc_randns = model.myosin.n + model.actin.n;
    std::vector<double> acc_rand(n_acc_randns);
    for (int i = 0; i < n_acc_randns; i++) {
        acc_rand[i] = gsl_rng_uniform(rng);
    }

    int offset = model.myosin.n * 3;
    double D = D_myosin;
    // Update myosin particles.
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        double dx = model.myosin.force[i].x * beta * D * dt +
                    model.myosin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[i * 3];
        double dy = model.myosin.force[i].y * beta * D * dt +
                    model.myosin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[i * 3 + 1];
        double dtheta = model.myosin.angular_force[i] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[i * 3 + 2] * M_PI / 5;
        model.myosin.displace(i, dx, dy, dtheta);
    }
    D = D_actin;
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        double dx = model.actin.force[i].x * beta * D * dt +
                    model.actin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 3];
        double dy = model.actin.force[i].y * beta * D * dt +
                    model.actin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 3 + 1];
        double dtheta = model.actin.angular_force[i] * beta * D * dt +
                        sqrt(2 * D * dt) * noise[offset + i * 3 + 2] * M_PI / 5;

        // (Optional debug printing if the displacement is large.)
        if (dx > 0.04 || dy > 0.04) {
            printf("actin %d\n", i);
            printf("dx: %f dy: %f\n", dx, dy);
            printf("force: %f %f\n", model.actin.force[i].x, model.actin.force[i].y);
            printf("force*beta*D*dt: %f %f\n",
                   model.actin.force[i].x * beta * D * dt,
                   model.actin.force[i].y * beta * D * dt);
            printf("velocity*dt: %f %f\n",
                   model.actin.velocity[i].x * dt, model.actin.velocity[i].y * dt);
            printf("sqrt(2*D*dt)*noise: %f %f\n",
                   sqrt(2 * D * dt) * noise[offset + i * 3],
                   sqrt(2 * D * dt) * noise[offset + i * 3 + 1]);
            printf("noise: %f %f\n",
                   noise[offset + i * 3], noise[offset + i * 3 + 1]);
        }
        model.actin.displace(i, dx, dy, dtheta);
    }
}
