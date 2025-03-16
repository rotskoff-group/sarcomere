#include "langevin.h"
#include "sarcomere.h"
#include "components.h"
#include "cxxopts.hpp"
#include <vector>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

int main(int argc, char* argv[]){
    int nsteps;
    int seed;
    double dt;
    double beta;
    double actin_diff_coeff;
    double myosin_diff_coeff;
    int update_dt_every;
    int update_myosin_every;
    int save_every;
    double k_on;
    double k_off;
    double base_lifetime;
    double lifetime_coeff;
    double k_aa;
    double kappa_aa;
    double k_am;
    double kappa_am;
    double v_am;
    int n_actins;
    int n_myosins;
    double Lx;
    double Ly;
    double actin_length;
    double myosin_length;
    double myosin_radius;
    double crosslinker_length;
    bool resume;
    bool directional;
    int n_fixed_myosins;

    std::string filename;
    std::string init_struc;
    try {
        cxxopts::Options options("sarcomere", "simulate sarcomere assembly using Monte Carlo");

        options.add_options()
            ("nsteps", "Number of steps", cxxopts::value<int>(nsteps)->default_value("500000"))
            ("seed", "Seed value", cxxopts::value<int>(seed)->default_value("0"))
            ("dt", "Time step", cxxopts::value<double>(dt)->default_value("0.00001"))
            ("beta", "Beta value", cxxopts::value<double>(beta)->default_value("241.0"))
            ("actin_diff_coeff", "Actin diffusion coefficient", cxxopts::value<double>(actin_diff_coeff)->default_value("0.1"))
            ("myosin_diff_coeff", "Myosin diffusion coefficient", cxxopts::value<double>(myosin_diff_coeff)->default_value("0.01"))
            ("update_dt_every", "Update dt every", cxxopts::value<int>(update_dt_every)->default_value("500"))
            ("update_myosin_every", "Update myosin every", cxxopts::value<int>(update_myosin_every)->default_value("1"))
            ("save_every", "Save every", cxxopts::value<int>(save_every)->default_value("200"))
            ("k_on", "k_on", cxxopts::value<double>(k_on)->default_value("1000"))
            ("k_off", "k_off", cxxopts::value<double>(k_off)->default_value("1"))
            ("base_lifetime", "Base lifetime", cxxopts::value<double>(base_lifetime)->default_value("0.001"))
            ("lifetime_coeff", "Lifetime coefficient", cxxopts::value<double>(lifetime_coeff)->default_value("0.4"))
            ("k_aa", "k_aa", cxxopts::value<double>(k_aa)->default_value("300"))
            ("kappa_aa", "kappa_aa", cxxopts::value<double>(kappa_aa)->default_value("50"))
            ("k_am", "k_am", cxxopts::value<double>(k_am)->default_value("50"))
            ("kappa_am", "kappa_am", cxxopts::value<double>(kappa_am)->default_value("50"))
            ("v_am", "v_am", cxxopts::value<double>(v_am)->default_value("5"))
            ("n_actins", "Number of actins", cxxopts::value<int>(n_actins)->default_value("50"))
            ("n_myosins", "Number of myosins", cxxopts::value<int>(n_myosins)->default_value("4"))
            ("Lx", "Lx", cxxopts::value<double>(Lx)->default_value("10"))
            ("Ly", "Ly", cxxopts::value<double>(Ly)->default_value("10"))
            ("actin_length", "Actin length", cxxopts::value<double>(actin_length)->default_value("1"))
            ("myosin_length", "Myosin length", cxxopts::value<double>(myosin_length)->default_value("1.5"))
            ("myosin_radius", "Myosin radius", cxxopts::value<double>(myosin_radius)->default_value("0.2"))
            ("crosslinker_length", "Crosslinker length", cxxopts::value<double>(crosslinker_length)->default_value("0.06"))
            ("resume", "Resume", cxxopts::value<bool>(resume)->default_value("false"))
            ("directional", "Directional", cxxopts::value<bool>(directional)->default_value("true"))
            ("n_fixed_myosins", "Number of fixed myosins", cxxopts::value<int>(n_fixed_myosins)->default_value("0"))
            ("filename", "Filename", cxxopts::value<std::string>(filename)->default_value("data/traj.h5"))
            ("initial_structure", "Type of initial structure", cxxopts::value<std::string>(init_struc)->default_value("random"))
            ("h, help", "Print usage");

        auto result = options.parse(argc, argv);

        if (result.count("help")) {
            std::cout << options.help() << std::endl;
            return 0;
        }
    } catch (const std::exception& e)
	  {
		    std::cerr << "Error parsing options: " << e.what() << std::endl;
		    return 1;
	  }


    gsl_rng * rng = gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set(rng,seed);
    std::vector <double> box(2);
    box[0] = Lx;
    box[1] = Ly;
    double diff_coeff_ratio = actin_diff_coeff/myosin_diff_coeff;
    Sarcomere model(n_actins, n_myosins, box, actin_length, myosin_length,
                        myosin_radius, crosslinker_length, k_on, k_off,
                        base_lifetime, lifetime_coeff, diff_coeff_ratio,
                          k_aa, kappa_aa, k_am, kappa_am, v_am,
                        filename,rng, seed, n_fixed_myosins, dt, directional);
    Langevin sim(model, beta, dt, actin_diff_coeff, myosin_diff_coeff, update_myosin_every, update_dt_every, save_every, resume);
    if (!resume){
        if (init_struc == "sarcomere") {
        sim.model.sarcomeric_structure();}
        // else if (init_struc == "sarcomere"){
        //     sim.model.sarcomeric_structure();
        // }
        else if (init_struc == "n_fixed_myosins"){
            sim.model.myosin_on_a_lattice();
        }
        else if (init_struc == "partial"){
            sim.model.partial_fix(n_fixed_myosins);
        }
        else if (init_struc == "cb"){
            sim.model.cb();
        }
        else if (init_struc == "bad_cb"){
            sim.model.bad_cb();
        }
    }
    if (!resume){
        sim.volume_exclusion(5000, rng,n_fixed_myosins);
    }
    sim.run_langevin(nsteps, rng, n_fixed_myosins);
    gsl_rng_free(rng);
    return 0;
}
