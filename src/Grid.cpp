#include "Grid.h"
#include "Config.h"
#include "common_functions.h"
#include "Cell.h"
#include "Node.h"
#include "Shipment_handler.h"

#include <iostream>
#include <ctime>
#include <chrono>
#include <fstream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <set>


Grid::Grid(const Config& C) :
    n_cells_1d(C.n_cells_1d), C(C), LS(C), current_timestep(0)
{
    readNodes(C.node_file);
    if(C.shipments_on)
    {
        SH = new Shipment_handler(C.shipment_file, node_map);
    }

    if(C.cell_method == 1)
    {
        makeCellsFixed(); //Construct static grid
        simple_distances = true;
    }
    else if(C.cell_method == 2)
    {
        makeCellsDynamicSquare(C.cell_max_farms); //Construct adaptive grid.
        simple_distances = false;
    }
    else
    {
        std::cout << "Incorrect cell creation method (" << C.cell_method
                  << "). Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
    for(size_t i = 0; i < cells.size(); i++)
    {
        cells[i]->set_distance_index(i);
    }

    makeShortestDistances(); //Precalculate all shortest distances between cells.
    initCells(); //Set up the cells' nodes' statuses.
    makeKernelValues(); //Precalculate the distance kernel values.

    if(C.opt_method == 3 or C.opt_method == 4 or C.opt_method == 5)
    {
        makeBinomialLookup(); //Generate lookup tables for binomial distributions.
    }
    initSeed(); //Set up the chosen seed scheme.
    all_tra_nodes.reserve(all_nodes.size());
}

Grid::~Grid()
{
    for(Cell* c : cells)
    {
        delete c;
    }
    for(Node* n : all_nodes)
    {
        delete n;
    }
    if(C.shipments_on)
    {
        delete SH;
    }
}

//Read landscape file with nodes.
void Grid::readNodes(std::string node_file)
{
    std::cout << "Loading nodes..." << std::flush;
    TimePoint process_start = common::get_wall_time();
    std::vector<double> all_x, all_y; //To store x and y values for finding boundaries

    std::ifstream f(node_file, std::ifstream::in);
    if(!f)
    {
        std::cout << "Node file " << node_file << " not found. Exiting..."
                  << std::endl; exit(EXIT_FAILURE);
    }

    if(f.is_open())
    {
        common::skipBOM(f);
        unsigned int n_lines = common::get_n_lines(f);
        all_nodes.reserve(n_lines);
        all_x.reserve(n_lines);
        all_y.reserve(n_lines);
        std::set<unsigned int> unique_groups;

        std::string line;
        while(std::getline(f, line))
        {
            std::vector<std::string> line_vector = common::split(line, '\t');
            if(!line_vector.empty())
            {
                unsigned int id = common::stringToNum<unsigned int>(line_vector[0]);
                unsigned int group = common::stringToNum<unsigned int>(line_vector[1]);
                double x = common::stringToNum<double>(line_vector[2]);
                double y = common::stringToNum<double>(line_vector[3]);
                unsigned int num_sp1 = common::stringToNum<unsigned int>(line_vector[4]);
                unsigned int num_sp2 = common::stringToNum<unsigned int>(line_vector[5]);
                Node* n = new Node(x, y, num_sp1, num_sp2, id, group, C, this);
                all_nodes.push_back(n);
                node_map[id] = n;
                all_x.push_back(x);
                all_y.push_back(y);
                group_node_map[group].push_back(n);
                unique_groups.insert(group);
            }
        }
        findBoundaries(all_x, all_y);
        groups.assign(unique_groups.begin(), unique_groups.end());
    }

    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done. " << all_nodes.size() << " nodes read from " << node_file
              << " (" << process_duration << " ms)" << std::endl;

    if(C.local_kernel == 1 or C.local_kernel == 3 or C.local_kernel == 4 or
       C.local_kernel == 5 or C.local_kernel == 6 or C.local_kernel == 7 or
       C.local_kernel == 9)
    {
        initTransmissionParametersUSDOS2();
    }
    else if(C.local_kernel == 2)
    {
        initTransmissionParametersGaussian();
    }
    else if(C.local_kernel == 8)
    {
        initTransmissionParametersHayama();
    }
    else
    {
        std::cout << "No initiation for kernel " << C.local_kernel << "? Please handle this "
                  << "explicitly around line " << __LINE__ << " in " << __FILE__ << "."
                  << std::endl; exit(EXIT_FAILURE);
    }
}

void Grid::initTransmissionParametersUSDOS2()
{
    std::cout << "Calculating transmission parameters..." << std::flush;
    TimePoint process_start = common::get_wall_time();

//    //Species 1
//    unsigned int sum_num1 = 0;
//    double sum_num_p1 = 0.0;
//    double sum_num_q1 = 0.0;
//    //Species 2
//    unsigned int sum_num2 = 0;
//    double sum_num_p2 = 0.0;
//    double sum_num_q2 = 0.0;
//
//    for(Node* node : all_nodes)
//    {
//        //Species 1
//        sum_num1 += node->get_num_sp1();
//        sum_num_q1 += std::pow(node->get_num_sp1(), C.q1);
//        sum_num_p1 += std::pow(node->get_num_sp1(), C.p1);
//        //Species 2
//        sum_num2 += node->get_num_sp2();
//        sum_num_q2 += std::pow(node->get_num_sp2(), C.q2);
//        sum_num_p2 += std::pow(node->get_num_sp2(), C.p2);
//    }
//
//    //Species 1
//    if(sum_num1 > 0) //If there arent any animals of this type the divisor will be 0 and the result nan.
//    {
//        norm_sus1 = sum_num1 / sum_num_q1;
//        norm_inf1 = sum_num1 / sum_num_p1;
//    }
//
//    //Species 2
//    if(sum_num2 > 0) //If there arent any animals of this type the divisor will be 0 and the result nan.
//    {
//        norm_sus2 = sum_num2 / sum_num_q2;
//        norm_inf2 = sum_num2 / sum_num_p2;
//    }

    for(Node* node : all_nodes)
    {
        node->set_susceptibility(C.S1 * std::pow(node->get_num_sp1(), C.q1) +
                                 C.S2 * std::pow(node->get_num_sp2(), C.q2));
        node->set_infectiousness(C.T1 * std::pow(node->get_num_sp1(), C.p1) +
                                 C.T2 * std::pow(node->get_num_sp2(), C.p2));
    }


    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done " << "(" << process_duration << " ms)." << std::endl;
}

void Grid::initTransmissionParametersGaussian()
{
    for(Node* node : all_nodes)
    {
        node->set_susceptibility(1.0);
        node->set_infectiousness(1.0);
    }
}

void Grid::initTransmissionParametersHayama()
{
    for(Node* node : all_nodes)
    {
        node->set_susceptibility(std::log(node->get_num_sp1() + node->get_num_sp2()));
        node->set_infectiousness(10.0 * std::log(node->get_num_sp1() + node->get_num_sp2()));
    }
}

void Grid::reset()
{
    seeded_nodes.clear();
    replicate_kc = 0;
    max_n_inf = 0;
    current_timestep = 0;
    for(Cell* c : cells)
    {
        c->reset();
    }
    group_inf_node_map.clear();
}

//Seed new nodes at the start of a replicate.
std::vector<Node*> Grid::seed()
{
    static std::vector<Node*>::iterator it = seed_vector.begin();
    seeded_nodes.clear();
    seeded_nodes.reserve(C.seed_n);
    for(size_t i = 0; i<C.seed_n; i++)
    {
        Node* seeded_node = *it;
        seeded_node->become_seeded();
        seeded_nodes.push_back(seeded_node);
        it++;
    }
    std::cout << seeded_nodes.size() << " nodes seeded." << std::endl;
    return seeded_nodes;
}

std::vector<Node*> Grid::get_seeded_nodes()
{
    return seeded_nodes;
}

//Move time one step forward and perform transmission calculations.
unsigned long long int Grid::incrementTime()
{
    current_timestep += 1;
    all_tra_nodes.clear();
    for(Cell* c : cells)
    {
        std::vector<Node*> cell_tra = c->get_transmitting_nodes();
        all_tra_nodes.insert(all_tra_nodes.end(), cell_tra.begin(), cell_tra.end());
    }

    //Do infections
    unsigned long long int timestep_kc = 0;
    if(C.opt_on_node_level)
    {
        timestep_kc = this->makeTransmissionNode();
    }
    else
    {
        timestep_kc = this->makeTransmissionCell();
    }
    replicate_kc += timestep_kc;
    if(C.shipments_on)
    {
        SH->make_transmission(current_timestep);
    }
    //Update all nodes' statuses.
    for(Cell* c : cells)
    {
        c->commitTimeStep();
    }

    std::vector<size_t> current_status = this->get_current_status();
    size_t current_n_inf = current_status[1] + current_status[2] + current_status[3];
    if(current_n_inf > max_n_inf)
        max_n_inf = current_n_inf;

    return timestep_kc; //Return the number of kernel requests.
}


void Grid::add_to_group_inf(Node* n, unsigned int group)
{
    group_inf_node_map[group].push_back(n);
}

std::vector<size_t> Grid::get_current_status()
{
    //Update all nodes and get summary data.
    std::vector<size_t> current_status_v = {0, 0, 0, 0, 0};
    for(Cell* c : cells)
        c->current_status(current_status_v);

    return current_status_v;
}

std::unordered_map<unsigned int, std::vector<Node*>> Grid::get_inf_nodes_by_group()
{
    return group_inf_node_map;
}


size_t Grid::get_time()
{
    return current_timestep;
}

void Grid::findBoundaries(std::vector<double>& all_x, std::vector<double>& all_y)
{
    //Setup variables to be used
    x_bounds.resize(2);
    y_bounds.resize(2);
    x_bounds[0] = all_x[0];
    x_bounds[1] = all_x[0];
    y_bounds[0] = all_y[0];
    y_bounds[1] = all_y[0];

    //Find smallest and largest x
    for(double x : all_x)
    {
        if(x <= x_bounds[0])
            x_bounds[0] = x;
        if(x >= x_bounds[1])
            x_bounds[1] = x;
    }

    //...and y
    for(double y : all_y)
    {
        if(y <= y_bounds[0])
            y_bounds[0] = y;
        if(y >= y_bounds[1])
            y_bounds[1] = y;
    }

    //Add a small amount to xmax and ymax so no nodes end up on the edge.
    x_bounds[1] += std::abs(x_bounds[1]) * 0.00001;
    y_bounds[1] += std::abs(y_bounds[1]) * 0.00001;

    //Get side length from boundaries
    double x_len = x_bounds[1] - x_bounds[0];
    double y_len = y_bounds[1] - y_bounds[0];
    if(x_len > y_len)
        side = x_len;
    else
        side = y_len;
    longest_distance = std::sqrt(2*(side*side));
    std::cout << "Grid boundaries set to x:[" << x_bounds[0] << ", " << x_bounds[1]
                 << "] and y:[" << y_bounds[0] << ", " << y_bounds[1] << "]. Side set to "
                 << side << ". Longest possible distance set to " << longest_distance << "."
                 << std::endl;
}

unsigned int Grid::get_n_cells()
{
    return (unsigned int)cells.size();
}

unsigned int Grid::get_n_groups()
{
    return group_node_map.size();
}

std::vector<unsigned int> Grid::get_groups()
{
    return groups;
}

size_t Grid::get_max_n_inf()
{
    return max_n_inf;
}

std::vector<Node*> Grid::get_all_nodes()
{
    return all_nodes;
}

double Grid::get_cell_size_mean()
{
    std::vector<unsigned int> sizes;
    sizes.reserve(cells.size());
    for(Cell* c : cells)
    {
        sizes.push_back(c->get_n_nodes());
    }
    return common::average(sizes);
}

double Grid::get_cell_size_variance_from(double val)
{
    std::vector<unsigned int> sizes;
    sizes.reserve(cells.size());
    for(Cell* c : cells)
    {
        sizes.push_back(c->get_n_nodes());
    }
    return common::variance_from(sizes, val);
}

double Grid::get_cell_size_stdev_from(double val)
{
    return std::sqrt(get_cell_size_variance_from(val));
}

size_t Grid::get_cell_cell_distance(size_t c1_d_idx, size_t c2_d_idx)
{
    return cellCellDistanceMatrix.at(c1_d_idx).at(c2_d_idx);
}

//Over-estimated infection probability between two cells.
double Grid::get_cell_cell_pover(size_t c1_d_idx, size_t c2_d_idx)
{
    return cellCellPoverMatrix.at(c1_d_idx).at(c2_d_idx);
}

//Kernel calls for the entire replicate.
unsigned long long int Grid::get_replicate_kc()
{
    return replicate_kc;
}

//Calls python code externally. Only if active in config file.
void Grid::plotGridAndNodes(bool draw_empty_cells)
{
    //create out file
    std::cout << "Creating a colorful plot of the grid...";
    bool grid_f_success = false;
    std::ofstream ofs_grid("_grid.dat", std::ofstream::out);
    if(ofs_grid.is_open())
    {
        grid_f_success = true;
        //get corners of each cell
        std::ostringstream oss_grid;
        for(Cell* c : cells)
        {
            if(draw_empty_cells == true or
               c->get_n_nodes() > 0)
            {
                std::vector<Point> corners = c->get_corners();
                //write to sstream
                for(Point p : corners)
                {
                    oss_grid << p.get_x();
                    oss_grid << "\t";
                    oss_grid << p.get_y();
                    oss_grid << "\t";
                }
                oss_grid << c->get_id();
                oss_grid << "\n";
            }
        }
        //write to file
        ofs_grid << oss_grid.str();
        //close file
        ofs_grid.close();
    }

    //nodes
    //create out file
    bool node_f_success = false;
    std::ofstream ofs_node("_nodes.dat", std::ofstream::out);
    if(ofs_node.is_open())
    {
        node_f_success = true;
        //get pos of each node
        std::ostringstream oss_node;
        for(Node* n : all_nodes)
        {
            Point position = n->get_position();
            //write to sstream
            oss_node << position.get_x();
            oss_node << "\t";
            oss_node << position.get_y();
            oss_node << "\t";
            oss_node << n->get_parent()->get_id();
            oss_node << "\n";
        }

        //write to file
        ofs_node << oss_node.str();
        //close file
        ofs_node.close();
    }

    //Plot using the .dat files
    system("python3 plot_grid_node_relation.py");

    //Delete the temporary files
    if(grid_f_success)
    {
        std::remove("_grid.dat");
    }
    if(node_f_success)
    {
        std::remove("_nodes.dat");
    }
    std::cout << "done." << std::endl;
}

//Calls python code externally. Only if active in config file.
void Grid::plotCellSizeHist(unsigned int max_nodes)
{
    //create out file
    std::cout << "Creating a histogram of cell sizes...";
    bool size_f_success = false;
    std::ofstream ofs_cells("_cellsize.dat", std::ofstream::out);
    if(ofs_cells.is_open())
    {
        size_f_success = true;
        //get corners of each cell
        std::ostringstream oss_cells;
        oss_cells << max_nodes << "\n";
        for(Cell* c : cells)
        {
            oss_cells << c->get_n_nodes() << "\n";
        }

        //write to file
        ofs_cells << oss_cells.str();
        //close file
        ofs_cells.close();
    }

    //Plot using the .dat files
    system("python3 plot_cell_hist.py");

    //Delete the temporary files
    if(size_f_success)
    {
        std::remove("_grid.dat");
    }
    std::cout << "done." << std::endl;
}

//Create a static grid.
void Grid::makeCellsFixed()
{
    std::cout << "Constructing fixed-size cells..." << std::flush;
    TimePoint process_start = common::get_wall_time();
    double cell_side = side / n_cells_1d;
    std::vector<Node*> remaining_nodes = all_nodes;
    std::vector<Node*> temp_node_v;
    temp_node_v.reserve(all_nodes.size());
    int added_nodes = 0;

    for(unsigned int i = 0; i < n_cells_1d; i++)
    {
        for(unsigned int j = 0; j < n_cells_1d; j++)
        {
            double ll_x = x_bounds[0] + j*cell_side;
            double ll_y = y_bounds[0] + i*cell_side;
            Cell* c = new Cell(ll_x, ll_y, cell_side, this, C);
            cells.push_back(c);

            //Add the nodes that belong to this cell.
            for(Node* n : remaining_nodes)
            {
                if(c->isWithin(n))
                {
                    c->addNode(n);
                    added_nodes++;
                }
                else
                {
                    temp_node_v.push_back(n);
                }

            }
            remaining_nodes.swap(temp_node_v); //Save the remaining nodes in rem_nodes vector.
            temp_node_v.clear(); //Empty the temp vector for next iteration.
        }
    }
    if(remaining_nodes.size() > 0)
    {
        std::cout << "All nodes were not added to cells during fixed cell creation. Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }

    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done. " << cells.size() << " identical cells created. (" << process_duration <<
                 " ms)" << std::endl;
}

//Make a adaptive grid by recursively dividing one large cell into 4 children until limit (lambda) is reached.
void Grid::makeCellsDynamicSquare(unsigned int max_nodes)
{
    std::cout << "Constructing dynamically sized square cells with max " << max_nodes
              << " nodes in each..." << std::flush;
    TimePoint process_start = common::get_wall_time();

    Cell supercell(x_bounds[0], y_bounds[0], side, this, C);
    supercell.addNodes(all_nodes);
    cells = supercell.subdivideQuad(max_nodes, true);

    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done. " << cells.size() << " cells created. (" << process_duration <<
                 " ms)" << std::endl;
}

//Precalculate distances between all grid cells.
void Grid::makeShortestDistances()
{
    std::string method_str;
    if(simple_distances)
        method_str = "simple";
    else
        method_str = "complex";

    std::cout << "Precalculating integer distances (rounded up) between all " <<  cells.size()
              << " cells using " << method_str << " function..." << std::flush;
    TimePoint process_start = common::get_wall_time();

    //Initiate the distance matrix.
    cellCellDistanceMatrix.resize(cells.size());
    for(size_t i = 0; i < cellCellDistanceMatrix.size(); i++)
    {
        cellCellDistanceMatrix.at(i).resize(cells.size(), 0);
    }

    //Find every cell-to-cell distance and save as an integer.
    for(Cell* c1 : cells)
    {
        for(Cell* c2 : cells)
        {
            double d = 0.0;
            if(simple_distances)
                d = c1->calc_shortestDistanceSimple(c2);
            else
                d = c1->calc_shortestDistanceComplex(c2);

            cellCellDistanceMatrix[c1->get_distance_index()][c2->get_distance_index()] = (unsigned int)(d + 0.5);
        }
    }

    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done. (" << process_duration << " ms)" << std::endl;
}

//Precalculate distance kernel values.
void Grid::makeKernelValues()
{
    std::cout << "Precalculating kernel values for all integer distances between 0 and "
              << int(longest_distance + 0.5) << "..." << std::flush;
    TimePoint process_precalc_kernel = common::get_wall_time();

    LS.init_kernel_lookup(longest_distance);

    unsigned long long int process_precalc_kernel_duration = common::delta_duration(process_precalc_kernel);
    std::cout << "done (" << process_precalc_kernel_duration << " ms)." << std::endl;

    std::cout << "Precalculating overestimated probability of infection between all " <<  cells.size()
              << " cells..." << std::flush;
    TimePoint process_pover_start = common::get_wall_time();

    //Initiate the p-over matrix.
    cellCellPoverMatrix.resize(cells.size());
    for(size_t i = 0; i < cellCellPoverMatrix.size(); i++)
    {
        cellCellPoverMatrix.at(i).resize(cells.size(), 0);
    }

    for(Cell* inf_cell : cells)
    {
        for(Cell* sus_cell : cells)
        {
            double p_over = LS.calc_infection_p_over(inf_cell, sus_cell);
            cellCellPoverMatrix.at(inf_cell->get_distance_index()).at(sus_cell->get_distance_index()) = p_over;
        }
    }
    unsigned long long int process_pover_duration = common::delta_duration(process_pover_start);
    std::cout << "done. (" << process_pover_duration << " ms)" << std::endl;
}

void Grid::makeBinomialLookup()
{
     std::cout << "Setting up binomial random variables lookup table " << std::flush;
    TimePoint process_start = common::get_wall_time();
    //Find number of nodes in largest cell.
    size_t largest_n = 1;
    for(Cell* c : cells)
    {
        size_t n = c->get_n_nodes();
        if(n > largest_n)
            largest_n = n;
    }
    std::cout << "(largest n: " << largest_n << ")..." << std::flush;
    LS.init_binomial_lookup(largest_n);
    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "done (" << process_duration << " ms)." << std::endl;
}


void Grid::initCells()
{
    for(Cell* c : cells)
    {
        c->init_node_statuses();
        c->init_transmission_parameters();
    }
}

void Grid::initSeed()
{
    if(C.seed_method == 1)
    {
        common::random_unique(all_nodes, C.n_replicates*C.seed_n, seed_vector);
    }

    else if(C.seed_method == 2)
    {
//        std::cout << "Seeding from file is currently unavailable. "
//                  << "Function unfinished." << std::endl;
//        ping();
//        exit(EXIT_FAILURE);
        //Read seed info from file and save to seed_vector.
        std::ifstream f(C.seed_file, std::ifstream::in);
        if(!f)
        {
            std::cout << "Seed file not found. Exiting..."
                      << std::endl; exit(EXIT_FAILURE);
        }
        if(f.is_open())
        {
            common::skipBOM(f);
            unsigned int n_lines = common::get_n_lines(f);
            std::vector<Node*> seed_vector_one_rep;
            std::string line;
            size_t line_n = 1;
            while(std::getline(f, line))
            {
                size_t seed_id = std::stoi(line);
                Node* seed_node;
                try
                {
                    seed_node = node_map.at(seed_id);
                }
                catch(const std::out_of_range& oorange_exception)
                {
                    std::cerr << "Error: " << oorange_exception.what() << std::endl
                              << seed_id << " in the " << C.seed_file
                              << " not found among farms." << std::endl;
                    exit(EXIT_FAILURE);
                }
                seed_vector_one_rep.push_back(seed_node);
                line_n++;
            }
            f.close();
            seed_vector.reserve(C.n_replicates * n_lines);
            for(size_t i = 0; i < C.n_replicates; i++)
            {
                seed_vector.insert(seed_vector.end(), seed_vector_one_rep.begin(), seed_vector_one_rep.end());
            }

        }
    }
    else if(C.seed_method == 3)
    {
        seed_vector.assign(C.n_replicates, node_map.at(C.specific_seed_id)); //This option is reused as id of node to seed in case of single spec. node seeding.
    }
    else if(C.seed_method == 4)
    {
        seed_vector.reserve(this->get_n_groups());
        //Seeds one random node in each group. The number of replicates is changed in main.cpp to equal n_groups
        for(auto& group_nodes_pair : group_node_map)
        {
            int idx = common::rand_int(0, group_nodes_pair.second.size() - 1);
            seed_vector.push_back(group_nodes_pair.second.at(idx));
        }
    }
    else
    {
        std::cout << "Unknown seed method (" << C.seed_method << "). Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }
}

//Performs transmission between nodes and cells.
unsigned long long int Grid::makeTransmissionNode()
{
    unsigned long long int timestep_kc = 0; //Keeps track of total number of kernel calls made this timestep
    unsigned long long int node_level_kc;


    for(size_t node_i = 0; node_i < all_tra_nodes.size(); node_i++)
    {
        node_level_kc = 0;
        Node* inf_n = all_tra_nodes[node_i];
        for(Cell* c_sus : cells)
        {
            if(c_sus->get_n_sus() > 0)
            {
                node_level_kc += LS.makeInfections(inf_n, c_sus);
            }
        }
        timestep_kc = timestep_kc + node_level_kc;
    }
    return timestep_kc;
}

//Performs transmission between cells.
unsigned long long int Grid::makeTransmissionCell()
{
    unsigned long long int timestep_kc = 0; //Keeps track of total number of kernel calls made this timestep
    unsigned long long int cell_level_kc;


    for(Cell* c1 : cells)
    {
        if(c1->get_n_inf() > 0)
        {
            cell_level_kc = 0;
            for(Cell* c2 : cells)
            {
                cell_level_kc += LS.makeInfections(c1, c2);
            }
            timestep_kc = timestep_kc + cell_level_kc;
        }
    }
    return timestep_kc;
}
