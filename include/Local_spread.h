#ifndef LOCAL_SPREAD_H
#define LOCAL_SPREAD_H

#include <vector>
#include <unordered_map>
#include <random>
#include "Point.h"

class Config;
class Node;
class Cell;

/**
This class handles the calculations related to the local spread of infection.
*/
class Local_spread
{
    //typedef double (*k_fun_ptr)(double); ///Typedef for a pointer to a kernel function.
    //typedef double (*d_fun_ptr)(Point, Point); ///Typedef for a pointer to a distance function.
    typedef unsigned long long int (Local_spread::*o_fun_node_ptr)(Node*, Cell*); ///Typedef of a pointer to optimization method based on node to cell transmission.
    typedef unsigned long long int (Local_spread::*o_fun_cell_ptr)(Cell*, Cell*); ///Typedef of a pointer to optimization method based on cell to cell transmission.
    typedef double (Local_spread::*p_fun_ptr)(Node*, Node*); ///Typedef of a pointer to probability function for true p.
    typedef double (Local_spread::*po_cell_fun_ptr)(Cell*, Cell*); ///Typedef of a pointer to probability function for overestimated p (to enter a cell).

    public:
        Local_spread(const Config& C);
        ~Local_spread();
        void reset();
        void init_kernel_lookup(double longest_distance);
        void init_binomial_lookup(size_t largest_n);
        double callKernel(double d);
        unsigned long long int makeInfections(Node* inf, Cell* recipient_cell);
        unsigned long long int makeInfections(Cell* inf_cell, Cell* recipient_cell);
        double calc_infection_p(Node* inf_node, Node* sus_node);
        double calc_infection_p_over(Cell* inf_cell, Cell* sus_cell);

    private:
        const Config& C;
        o_fun_node_ptr o_function_node; ///Pointer to optimization function in use for node to cell transmission.
        o_fun_cell_ptr o_function_cell; ///Pointer to optimization function in use for cell to cell transmission.
        p_fun_ptr p_function; ///Pointer to probability function in use (depends on kernel choice).
        po_cell_fun_ptr po_cell_function; ///Pointer to probability function in use (for overestimated cell to cell probability).
        std::vector<double> kernel_lookup;
        std::vector<std::vector<std::binomial_distribution<size_t>*>*> binomial_lookup;
        std::mt19937* generator;
        std::vector<double> binomial_lookup_ps = {1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 1.0e-2, 1.0e-3,
                                                  1.0e-4, 1.0e-5, 1.0e-6, 1.0e-7, 1.0e-8, 5.0e-9};

        double generate_bin_rv_lookup(size_t n, double p, size_t& res);
        double calc_infection_p_USDOS2(Node* inf_node, Node* sus_node); ///Get a probability based on distance between and norm. sus and inf of sus_node and inf_node.
        double calc_infection_p_over_USDOS2(Node* inf_node, Cell* sus_cell);///Get p_over for USDOS2 binomial & keeling methods between inf node and cell.
        double calc_infection_p_over_USDOS2(Cell* inf_cell, Cell* sus_cell);///Get p_over for USDOS2 binomial & keeling methods between inf cell and cell.
        double calc_infection_p_USDOS2_lookup(Node* inf_node, Node* sus_node); ///Get a probability based on distance between and norm. sus and inf of sus_node and inf_node using a lookup table for kernel values.

        unsigned long long int localPairwise(Node* inf, Cell* recipient_cell); ///Performs transmission using PW alg. from infectious node to cell.
        unsigned long long int localPairwise(Cell* inf_cell, Cell* recipient_cell); ///Performs transmission using PW alg. from cell to cell.
        unsigned long long int localKeeling(Node* inf, Cell* recipient_cell); ///Performs transmission using CE alg. from infectious node to cell.
        unsigned long long int localBinomial(Node* inf, Cell* recipient_cell); ///Performs transmission using simple CS alg. from infectious node to cell.
        unsigned long long int localBinomial_CtoC(Cell* inf_cell, Cell* recipient_cell); ///Performs transmission using CS alg. from cell to cell.
};

#endif // LOCAL_SPREAD_H
