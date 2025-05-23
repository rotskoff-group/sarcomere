#ifndef SARCOMERE_H
#define SARCOMERE_H

#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <unordered_map>
#include <utility>
#include <omp.h>

#include "components.h"
#include "utils.h"
#include "interaction.h"
#include "geometry.h"
#include "h5_utils.h"
#include "neighborlist.h"

using vector = std::vector<double>;
using interaction = geometry::am_interaction;
const double EPS = 1e-6;

class Sarcomere {
    public:
        Filament actin;
        Myosin myosin;
        NeighborList neighbor_list;
        utils::MoleculeConnection myosinIndicesPerActin;
        utils::MoleculeConnection actinIndicesPerMyosin;
        utils::MoleculeConnection actinIndicesPerActin;

        std::vector<std::vector<int>> actin_actin_bonds, actin_actin_bonds_prev;
        vector box;
        double k_aa, kappa_aa, cb_mult_factor, k_on, k_off,
               kappa_am, k_am, v_am, crosslinker_length, skin_distance, total_energy, cutoff_radius;

        std::vector<std::vector<interaction>> am_interaction;
        vector actin_crosslink_ratio;
        std::vector<int> actin_n_bonds;
        std::vector<std::pair<std::vector<int>, std::vector<int>>> actin_neighbors_by_species;
        vector actin_basic_tension;
        
        gsl_rng* rng;
        std::string filename;

        std::vector<std::vector<vec>> actin_forces_temp;
        std::vector<std::vector<vec>> myosin_forces_temp, myosin_velocities_temp;
        std::vector<std::vector<double>> actin_angular_forces_temp, myosin_angular_forces_temp, actin_cb_strengths_temp;
        std::vector<utils::MoleculeConnection> actinIndicesPerMyosin_temp;
        std::vector<gsl_rng*> rng_engines;

        // Constructor & Destructor
        Sarcomere();
        Sarcomere(int& n_actins, int& n_myosins, vector box0, double& actin_length, double& myosin_length,
                  double& myosin_radius, double& crosslinker_length, double& k_on, double& k_off,
                  double& cb_mult_factor, double& k_aa, double& kappa_aa, double& k_am, double& kappa_am, double& v_am,
                  std::string& filename, gsl_rng* rng, int& seed);
        ~Sarcomere();

        // Methods
        void myosin_on_a_lattice();
        void partial_fix();
        void cb();
        void bad_cb();
        void sarcomeric_structure();
        void update_system();
        void new_file();
        void save_state();
        void load_state();

    private:
        void _update_neighbors();
        void _set_to_zero();
        void _process_actin_myosin_binding(int& i);
        void _process_catch_bonds(int& i);
        void _actin_myosin_force(int& i);
        void _volume_exclusion();
        void _myosin_repulsion(int& i, int& j);
        void _actin_repulsion(int& i, int& j);
        double _get_cb_strength(int& i, int& j);
        void _set_cb(int& i, int& j, double& normalized_strength, bool& add_connection);
        void _set_cb(int& i, std::vector<int> indices, vector cb_strength);
};

#endif // SARCOMERE_H
