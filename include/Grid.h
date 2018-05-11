#ifndef GRID_H
#define GRID_H

#include "Local_spread.h"

#include <vector>
#include <string>
#include <unordered_map>
#include <list>

class Config;
class Node;
class Cell;
class Shipment_handler;

/**
The grid class represents the whole area where the nodes and grid cells
are held.
*/
class Grid
{
    public:
        /**
        Construct using the side length of the entire grid and desired
        number of cells along that side as arguments.
        */
        Grid(const Config& C);
        ~Grid();

        void reset(); ///Resets all nodes to susceptible status.
        std::vector<Node*> seed(); ///Makes the next n number of nodes in seed_vector infectious
        std::vector<Node*> get_seeded_nodes(); ///Returns the nodes from the last seeding event.
        unsigned long long int incrementTime(); ///Moves time forward one timestep. Returns number of kernel calls required during that timestep. Uses openmp.
        void add_to_group_inf(Node* n, unsigned int group); ///Call to add a node to the map of infected nodes by group.
        std::vector<size_t> get_current_status();
        std::unordered_map<unsigned int, std::vector<Node*>> get_inf_nodes_by_group();
        size_t get_time(); ///Returns the current timestep.
        double get_cell_size_mean();
        double get_cell_size_variance_from(double val);
        double get_cell_size_stdev_from(double val);
        unsigned long long int get_replicate_kc();
        unsigned int get_n_cells(); ///Returns the current number of cells.
        unsigned int get_n_groups(); ///Returns the number of groups found in the node data.
        std::vector<unsigned int> get_groups(); ///Returns a vector of all unique groups in the node data.
        size_t get_max_n_inf();
        std::vector<Node*> get_all_nodes(); //Returns the vector of all nodes.
        size_t get_cell_cell_distance(size_t c1_d_idx, size_t c2_d_idx); ///Looks up the distance between the cells with the distance indices of the argumetns and returns it.
        double get_cell_cell_pover(size_t c1_d_idx, size_t c2_d_idx); ///Looks up and returns precalculated overestimated probability of infection between the two cells.
        void plotGridAndNodes(bool draw_empty_cells);
        void plotCellSizeHist(unsigned int max_nodes);
        void plotGridCellsAndDistances();
        void output();

    private:
        double side; ///Side length of entire grid.
        double longest_distance; ///Distance between two opposite corners or the grid.
        unsigned int n_cells_1d; ///Number of cells in one dimension.
        unsigned int max_n_inf = 0;
        const Config& C;
        Local_spread LS;
        Shipment_handler* SH;
        bool simple_distances = false;
        size_t current_timestep;
        double norm_inf1 = 1.0; ///Normalization factor for infectiousness, species 1.
        double norm_sus1 = 1.0; ///Normalization factor for susceptibility, species 1.
        double norm_inf2 = 1.0; ///Normalization factor for infectiousness, species 2.
        double norm_sus2 = 1.0; ///Normalization factor for susceptibility, species 2.
        std::vector<double> x_bounds, y_bounds; ///X and Y boundaries.
        std::vector<Node*> all_nodes; ///All nodes within the grid.
        std::vector<Node*> all_tra_nodes; ///The nodes that are currently transmitting infection. Updated every timestep in incrementTime().
        std::vector<Cell*> cells; ///All cells within the grid.
        std::vector<unsigned int> groups; ///All different unique groups found in the node input file.
        std::unordered_map<unsigned int, std::vector<Node*>> group_node_map; ///Gives access to all nodes belonging to each group by group id.
        std::unordered_map<unsigned int, std::vector<Node*>> group_inf_node_map; ///Gives access to all nodes that have become infected during the replicate belonging to each group by group id.
        std::vector<std::vector<unsigned int>> cellCellDistanceMatrix; ///Matrix containing all distances between cells rounded up to nearest int.
        std::vector<std::vector<double>> cellCellPoverMatrix; ///Matrix containing all precalculated overestimated probabilities of infection between cells.
        std::vector<Node*> seed_vector; ///Filled with nodes to be seeded.
        std::vector<Node*> seeded_nodes; ///Stores the nodes that were seeded in the beginning of the current replicate.
        std::unordered_map<unsigned int, Node*> node_map; ///Gives access to nodes by id.
        unsigned long long int replicate_kc = 0;

        void readNodes(std::string node_file); ///Reads and creates nodes from an input file.
        void initTransmissionParametersUSDOS2(); ///Calculates norm_inf and norm_sus and sets all farms' susceptibility and infectiousness.
        void initTransmissionParametersGaussian();
        void initTransmissionParametersHayama();
        void findBoundaries(std::vector<double>& all_x, std::vector<double>& all_y);
        void makeCellsFixed(); ///Creates a set of identical cells based on members n_cells_1d and side.
        void makeCellsDynamicSquare(unsigned int max_nodes); ///Creates a set of square cells with dynamic size.
        void makeShortestDistances(); ///Calculates the shortest distances between all cells.
        void makeKernelValues(); ///Precalculates the kernel values for all shortest distances between cells.
        void makeBinomials(); ///Sets up all the binomial random variables to use for every other cell.
        void makeBinomialLookup();
        void initCells(); ///Calls init() on each cell.
        void initSeed(); ///Reads seed file if applicable and fills seed vector with seed nodes.

        unsigned long long int makeTransmissionNode();
        unsigned long long int makeTransmissionCell();
};

#endif // GRID_H
