#ifndef CONFIG_H
#define CONFIG_H

#include "common_functions.h"

#include <string>
#include <vector>
#include <thread>

class Config
{
    public:
        Config(std::string config_file);
        ~Config();

        unsigned int max_timesteps; ///Maximum number of time steps before exiting.
        unsigned int n_replicates; ///Number of replicates to run with same grid.
        unsigned int opt_method; ///Method for speeding up disease spread calculations (0 = none, pairwise comp., 1 = Keeling, 2 = binomial)
        bool opt_on_node_level; ///Does the chosen optimiation method work on node to cell or cell to cell level?
        unsigned int lookup_table; ///Use a precalculated lookup table for the kernel evaluated at all integer distances.
        unsigned int cell_method; ///Cell creation method (1 = fixed, 2 = dynamic).
        unsigned int n_cells_1d; ///Number of cells along one axis of the grid.
        unsigned int cell_max_farms; ///Max farms before subdividing a cell.
        unsigned int inc_time; ///Incubation time, between exposure and infectiousness.
        unsigned int det_time; ///Detection time, between infectiousness and detection.
        unsigned int rem_time; ///Removal time. Delay between detected and removed (immune, vaccinated, culled).
        unsigned int seed_method; ///Seed method to use (1 = random, 2 = from file, 3 = one specific node).
        unsigned int seed_n; ///Number of nodes to seed each replicate.
        unsigned int specific_seed_id = 99999999; ///This is the node that will get seeded if seed_metod = 3.
        unsigned int local_kernel; ///Function to use for distance dependent probabilities (1 = USDOS2).
        unsigned int distance_method = 1; ///Function to use for calculating distances (1 = normal).
        std::string batch_name; ///Name of this batch of simulations
        std::string node_file; ///Name of file containing node data.
        std::string seed_file; ///Name of file with seed information.
        std::string shipment_file; ///Name of file containing a list of shipments with shipment id, origin id, destination id, timestep of shipment and number of animals; all tab-separated. One shipment/line.
        std::vector<int> output_levels; ///Output "level". Gridsim will output replicate results when the cumulative number of infected reaches these levels. Will also exit after the last output.
        bool make_plots; ///Output plots of the landscape.
        bool shipments_on = false; ///Is spread by shipments turned on.
        k_fun_ptr k_function; ///Pointer to kernel function in use.
        d_fun_ptr d_function; ///Pointer to distance function in use.

        /* Local spread infectiousness and susceptibility parameters.
        q1 & q2 - SUSCEPTIBILITY exponents for species 1 & 2.
        p1 & p2 - INFECTIOUSNESS (TRANSMISSIBILITY) exponents for sp 1 & 2.
        S1 & S2 - SUSCEPTIBILITY scale parameters for sp 1 & 2.
        T1 & T2 - INFECTIOUSNESS (TRANSMISSIBILITY) scale parameters for sp 1 & 2. */
        double q1 = 0.42;
        double p1 = 0.41;
        double S1 = 5.7;
        double T1 = 0.00082;
        double q2 = 0.42;
        double p2 = 0.41;
        double S2 = 5.7;
        double T2 = 0.00082;

        std::string as_str(); ///Returns a string representation of all settings.
    private:
        std::vector<std::string> parameter_vector;
};

#endif // CONFIG_H
