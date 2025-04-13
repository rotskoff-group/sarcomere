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
Langevin::Langevin(Sarcomere& model0, double& beta0, double& dt0, double& D0_actin_trans,
                    double& D0_actin_rot, double& D0_myosin_trans, double& D0_myosin_rot, int& save_every0, bool& resume)
    : model(model0), beta(beta0), dt(dt0), D_actin_trans(D0_actin_trans),
        D_actin_rot(D0_actin_rot),D_myosin_trans(D0_myosin_trans), D_myosin_rot(D0_myosin_rot),
      save_every(save_every0)
{
    if (resume) {
        int n_frames;
        model.load_state(n_frames);
        start_step = (n_frames-1) * save_every + 1;
        printf("Resuming from file %s at step %d\n", model.filename.c_str(), start_step);
        // // Print myosin details.
        // for (int i = 0; i < model.myosin.n; i++) {
        //     printf("Myosin %d: %f %f\n", i, model.myosin.center[i].x, model.myosin.center[i].y);
        //     printf("Myosin endpoints: (%f %f), (%f %f)\n",
        //            model.myosin.left_end[i].x, model.myosin.left_end[i].y,
        //            model.myosin.right_end[i].x, model.myosin.right_end[i].y);
        // }
        // // Print actin details.
        // for (int i = 0; i < model.actin.n; i++) {
        //     printf("Actin %d: %f %f\n", i, model.actin.center[i].x, model.actin.center[i].y);
        //     printf("Actin endpoints: (%f %f), (%f %f)\n",
        //            model.actin.left_end[i].x, model.actin.left_end[i].y,
        //            model.actin.right_end[i].x, model.actin.right_end[i].y);
    } else {
        start_step = 0;
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
void Langevin::run_langevin(int& nsteps, gsl_rng* rng, int& fix_myosin) {
    double start, end;
    for (int i = start_step; i < nsteps; i++) {
        if (i % save_every == 0) {
            std::cout << "Step " << i << std::endl;
            // Optionally protect saving with the mutex:
            // std::lock_guard<std::mutex> lock(save_mutex);
            model.save_state();
            start = omp_get_wtime();
        }
        model.update_system();
        sample_step(dt, rng, fix_myosin);
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
        sample_step(dt, rng, fix_myosin);
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
void Langevin::sample_step(double& dt, gsl_rng* rng, int& fix_myosin) {
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
    double D = D_myosin_trans;
    double D_rot = D_myosin_rot;
    // Update myosin particles.
    for (int i = fix_myosin; i < model.myosin.n; i++) {
        double dx = model.myosin.force[i].x * beta * D * dt +
                    model.myosin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[i * 3];
        double dy = model.myosin.force[i].y * beta * D * dt +
                    model.myosin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[i * 3 + 1];
        double dtheta = model.myosin.angular_force[i] * beta * D_rot * dt +
                        sqrt(2 * D_rot * dt) * noise[i * 3 + 2];
        // (Optional debug printing if the displacement is large.)
        if (dx > 0.06 || dy > 0.06) {
            printf("myosin %d\n", i);
            printf("dx: %f dy: %f\n", dx, dy);
            printf("force: %f %f\n", model.myosin.force[i].x, model.myosin.force[i].y);
            printf("force*beta*D*dt: %f %f\n",
                   model.myosin.force[i].x * beta * D * dt,
                   model.myosin.force[i].y * beta * D * dt);
            printf("velocity*dt: %f %f\n",
                     model.myosin.velocity[i].x * dt, model.myosin.velocity[i].y * dt);  
        }
        model.myosin.displace(i, dx, dy, dtheta);
    }
    // Update actin particles.
    for (int i = 0; i < model.actin.n; i++) {
        if (model.actin.cb_strength[i] > 1e-3){
            D = D_myosin_trans;
            D_rot = D_myosin_rot;
        }
        else{
            D = D_actin_trans;
            D_rot = D_actin_rot;
        }
        double dx = model.actin.force[i].x * beta * D * dt +
                    model.actin.velocity[i].x * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 3];
        double dy = model.actin.force[i].y * beta * D * dt +
                    model.actin.velocity[i].y * dt +
                    sqrt(2 * D * dt) * noise[offset + i * 3 + 1];
        double dtheta = model.actin.angular_force[i] * beta * D_rot * dt +
                        sqrt(2 * D_rot * dt) * noise[offset + i * 3 + 2];

        // (Optional debug printing if the displacement is large.)
        if (dx > 0.06 || dy > 0.06) {
            printf("actin %d\n", i);
            printf("dx: %f dy: %f\n", dx, dy);
            printf("force: %f %f\n", model.actin.force[i].x, model.actin.force[i].y);
            printf("force*beta*D*dt: %f %f\n",
                   model.actin.force[i].x * beta * D * dt,
                   model.actin.force[i].y * beta * D * dt);
            printf("velocity*dt: %f %f\n",
                   model.actin.velocity[i].x * dt, model.actin.velocity[i].y * dt);
        }
        model.actin.displace(i, dx, dy, dtheta);
    }
}
