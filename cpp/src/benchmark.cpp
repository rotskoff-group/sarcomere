// bench_2d.cpp
#include <benchmark/benchmark.h>
#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include <vector>
#include <gsl/gsl_rng.h>
#include <omp.h>

// Fixture for benchmarking run_langevin in 2D
class RunLangevin2DBenchmark : public benchmark::Fixture {
public:
    void SetUp(const ::benchmark::State& state) override {
        // 1) set OMP threads = first argument of state
        omp_set_num_threads(state.range(0));

        // 2) simulation / benchmark params (nsteps small)
        nsteps = 10;
        seed = 0;
        dt = 1e-5;
        beta = 241.0;
        actin_diff_t = 1;      actin_diff_r = 2;
        myosin_diff_t = 0.05;  myosin_diff_r = 0.05;
        save_every = 200;
        k_on = 100;  k_off = 1;
        base_lifetime = 0.001;  lifetime_coeff = 0.4;
        k_aa = 300;  kappa_aa = 50;
        k_am = 50;   kappa_am = 50;   v_am = 5;
        n_actins = 1500;  n_myosins = 100;
        Lx = 12;  Ly = 6;
        // no Lz in 2D
        actin_length = 1;   myosin_length = 1.5;
        myosin_radius = 0.2;  myosin_radius_ratio = 0.75;
        crosslinker_length = 0.06;
        resume = false;  directional = true;
        n_fixed_myosins = 0;
        filename = "traj2d.h5";
        init_struc = "random";
        boundary_condition = "periodic";
        

        // build 2D box
        std::vector<double> box = {Lx, Ly};

        // RNG
        rng = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng, seed);

        // model & sim
        double diff_coeff_ratio = actin_diff_t / myosin_diff_t;
        model = new Sarcomere(
            n_actins, n_myosins, box,
            actin_length, myosin_length,
            myosin_radius, myosin_radius_ratio,
            crosslinker_length, k_on, k_off,
            base_lifetime, lifetime_coeff, diff_coeff_ratio,
            k_aa, kappa_aa, k_am, kappa_am, v_am,
            filename, rng, seed, n_fixed_myosins, dt, directional, boundary_condition
        );
        sim = new Langevin(
            *model, beta, dt,
            actin_diff_t, actin_diff_r,
            myosin_diff_t, myosin_diff_r,
            save_every, resume
        );

        // initial structure & volume‐exclusion
        if (!resume) {
          if (init_struc == "sarcomere")  sim->model.sarcomeric_structure();
          else if (init_struc == "partial") sim->model.partial_fix(n_fixed_myosins);
          else if (init_struc == "cb")      sim->model.cb();
        }
        sim->volume_exclusion(1, rng, n_fixed_myosins);
    }

    void TearDown(const ::benchmark::State&) override {
        gsl_rng_free(rng);
        delete sim;
        delete model;
    }

    // members…
    int    nsteps, seed, save_every;
    int    n_actins, n_myosins, n_fixed_myosins;
    double dt, beta;
    double actin_diff_t, actin_diff_r;
    double myosin_diff_t, myosin_diff_r;
    double k_on, k_off, base_lifetime, lifetime_coeff;
    double k_aa, kappa_aa, k_am, kappa_am, v_am;
    double Lx, Ly, actin_length, myosin_length, myosin_radius, myosin_radius_ratio;
    double crosslinker_length;
    bool   resume, directional;
    std::string filename, init_struc, boundary_condition;

    gsl_rng* rng;
    Sarcomere* model;
    Langevin*  sim;
};

BENCHMARK_DEFINE_F(RunLangevin2DBenchmark, RunLangevin2D)(benchmark::State& st) {
    for (auto _ : st) {
        sim->run_langevin(nsteps, rng, n_fixed_myosins);
    }
    st.SetItemsProcessed(st.iterations());
}

// register for a range of thread‐counts: 2,4,8,16,32,64,128,256
BENCHMARK_REGISTER_F(RunLangevin2DBenchmark, RunLangevin2D)
    ->Arg(2)
    ->Arg(4)
    ->Arg(8)
    ->Arg(16)
    ->Arg(32)
    ->Arg(64)
    ->Arg(128)
    ->Arg(256);

BENCHMARK_MAIN();
