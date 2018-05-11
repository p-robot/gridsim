#include "Local_spread.h"
#include "common_functions.h"
#include "Config.h"
#include "Cell.h"
#include "Node.h"

#include <random>
#include <algorithm>

Local_spread::Local_spread(const Config& C) : C(C)
{
    //Setup optimization function pointers based on option set in config file.
    if(C.opt_method == 1)
    {
        o_function_node = &Local_spread::localPairwise;
        o_function_cell = nullptr;
    }
    else if(C.opt_method == 2)
    {
        o_function_node = &Local_spread::localKeeling;
        o_function_cell = nullptr;
    }
    else if(C.opt_method == 3)
    {
        o_function_node = &Local_spread::localBinomial;
        o_function_cell = nullptr;
    }
    else if(C.opt_method == 4)
    {
        o_function_node = nullptr;
        o_function_cell = &Local_spread::localBinomial_CtoC;
    }
    else
    {
        std::cout << "Unknown optimization method selected (" << C.opt_method
                  << ")." << std::endl << "Exiting..." << std::endl;
        exit(EXIT_FAILURE);
    }

    //Set pointers to functions to use for calculating probabilities.
    if(C.local_kernel == 1 or C.local_kernel == 2 or
       C.local_kernel == 3 or C.local_kernel == 4 or
       C.local_kernel == 5 or C.local_kernel == 6 or
       C.local_kernel == 7 or C.local_kernel == 8 or
       C.local_kernel == 9)
    {
        po_cell_function = &Local_spread::calc_infection_p_over_USDOS2; //Overload for cell & cell

        if(C.lookup_table == 0)
        {
            p_function = &Local_spread::calc_infection_p_USDOS2;
        }
        else if(C.lookup_table == 1)
        {
            p_function = &Local_spread::calc_infection_p_USDOS2_lookup;
        }
        else
        {
            std::cout << "Unknown kernel lookup table option (6) selected (" << C.local_kernel
                      << ")." << std::endl << "Exiting..." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    generator = new std::mt19937(std::chrono::system_clock::now().time_since_epoch().count());
    std::sort(binomial_lookup_ps.begin(), binomial_lookup_ps.end());
}

Local_spread::~Local_spread()
{
    //dtor
}

void Local_spread::init_kernel_lookup(double longest_distance)
{
    //Precalculate kernel values for all possible integer distances within the landscape.
    kernel_lookup.assign(int(longest_distance + 2.5), -1.0);
    for(double d = 0.0; d < longest_distance + 1.0; d+=1.0)
    {
        double k = this->callKernel(d);
        kernel_lookup.at(int(d+0.5)) = k;
    }
}

//Initiates a set of binomial random value generators for all combinations of p and n in binomial_lookup_ps and [1 to largest cell] respectively.
void Local_spread::init_binomial_lookup(size_t largest_n)
{
    for(size_t n=0; n<largest_n+1; n++)
    {
//        std::vector<Bin5tbl*>* v = new std::vector<Bin5tbl*>;
        std::vector<std::binomial_distribution<size_t>*>* v = new std::vector<std::binomial_distribution<size_t>*>;
        v->assign(binomial_lookup_ps.size(), nullptr);
        for(size_t i = 0; i < binomial_lookup_ps.size(); i++)
        {
            double p = binomial_lookup_ps[i];
//            Bin5tbl* b = new Bin5tbl(n, p, true);
            std::binomial_distribution<size_t>* b = new std::binomial_distribution<size_t>(n, p);
            v->at(i) = b;
        }
        binomial_lookup.push_back(v);
    }
}

//Performs the transmission between one infectious node and one specific cell based on chosen optimization function (CE or PW).
unsigned long long int Local_spread::makeInfections(Node* inf, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 0;
    if(inf->get_parent() == recipient_cell) //If it is within-cell spread, use pairwise.
        n_kernel_calls = localPairwise(inf, recipient_cell);
    else //Otherwise use the chosen opt function.
        n_kernel_calls = (this->*o_function_node)(inf, recipient_cell);
    return n_kernel_calls;
}

//Performs the transmission between one cell with infectious nodes and one other cell based on chosen optimization function (CS).
unsigned long long int Local_spread::makeInfections(Cell* inf_cell, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 0;
    if(inf_cell == recipient_cell) //If it is within-cell spread, use pairwise.
        n_kernel_calls = localPairwise(inf_cell, recipient_cell);
    else //Otherwise use the chosen opt function.
        n_kernel_calls = (this->*o_function_cell)(inf_cell, recipient_cell);
    return n_kernel_calls;
}

//Returns the true probability of infection between two nodes.
double Local_spread::calc_infection_p(Node* inf_node, Node* sus_node)
{
    return (this->*p_function)(inf_node, sus_node);
}

//Returns the over-estimated cell-to-cell probability of infection
double Local_spread::calc_infection_p_over(Cell* inf_cell, Cell* sus_cell)
{
    return (this->*po_cell_function)(inf_cell, sus_cell);
}

//Generates a binomial random variate based on n and p and saves it to res. However p is rounded to nearest p in binomial_lookup_ps and that rounded value is returned.
double Local_spread::generate_bin_rv_lookup(size_t n, double p, size_t& res)
{
    //Round p up to nearest p in the lookup table.
    auto const it = std::lower_bound(binomial_lookup_ps.begin(), binomial_lookup_ps.end(), p);
    size_t idx = it - binomial_lookup_ps.begin();
    res = binomial_lookup[n]->operator[](idx)->operator()(*generator);
    return *it;
}

///This function takes care of the calculation of distance between nodes and
///evaluating the kernel at that distance. Then it uses the normalized susceptibility
///and infectiousness of the two nodes together with the kernel value to return the
///probability of sus_node to get infected by inf_node.
double Local_spread::calc_infection_p_USDOS2(Node* inf_node, Node* sus_node)
{
    double d = C.d_function(inf_node->get_position(), sus_node->get_position());
    double k = callKernel(d);
    double inf_j = inf_node->get_infectiousness();
    double sus_i = sus_node->get_susceptibility();
    double exponent =  inf_j * sus_i * k * -1;
    return common::oneMinusExp(exponent);
}

///This function takes calculates an overestimated probability of infection
///between an infected farm and another cell. It uses the largest susceptibility
///of the receiving cell and the infectiousness of the infected node together
///with the kernel value to return the overestimated probability of 'entering' the cell.
double Local_spread::calc_infection_p_over_USDOS2(Node* inf_node, Cell* sus_cell)
{
    double d = double(inf_node->get_parent()->get_shortestDistanceTo(sus_cell));
    double k = kernel_lookup.at(d);
    double inf_j = inf_node->get_infectiousness();
    double sus_i = sus_cell->get_max_susceptibility();
    double exponent = inf_j * sus_i * k * -1;
    return common::oneMinusExp(exponent);
}

///This function takes calculates an overestimated probability of infection
///between two given cells based on the shortest distance between them as well
///as their maximum respective susceptiblility and infectiousness.
double Local_spread::calc_infection_p_over_USDOS2(Cell* inf_cell, Cell* sus_cell)
{
    double d = double(inf_cell->get_shortestDistanceTo(sus_cell));
    double k = this->callKernel(d);
    double inf_j = inf_cell->get_max_infectiousness();
    double sus_i = sus_cell->get_max_susceptibility();
    double exponent = inf_j * sus_i * k * -1;
    return common::oneMinusExp(exponent);
}

///This function takes care of the calculation of distance between nodes and
///fetches the pre-evaluated kernel value for that distance. Then it uses the
///normalized susceptibility and infectiousness of the two nodes together with
///the kernel value to return the probability of sus_node to get infected by inf_node.
double Local_spread::calc_infection_p_USDOS2_lookup(Node* inf_node, Node* sus_node)
{
    int d = int(C.d_function(inf_node->get_position(), sus_node->get_position()) + 0.5);
    double k = kernel_lookup.at(d);
    double inf_j = inf_node->get_infectiousness();
    double sus_i = sus_node->get_susceptibility();
    double exponent =  inf_j * sus_i * k * -1;
    return common::oneMinusExp(exponent);
}

double Local_spread::callKernel(double d)
{
    return C.k_function(d);
}

//This is the core of the PW algorithm.
unsigned long long int Local_spread::localPairwise(Node* inf, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 0;
    for(Node* sus : recipient_cell->get_nodes(0))
    {
        double p = (this->*p_function)(inf, sus);
        n_kernel_calls++;
        if(common::unif_rand() < p)
            sus->become_exposed();
    }
    return n_kernel_calls;
}

//PW for two entire cells at once.
unsigned long long int Local_spread::localPairwise(Cell* inf_cell, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 0;
    for(Node* inf : inf_cell->get_nodes(2)) //For every infectious node in cell 1.
    {
        n_kernel_calls += this->localPairwise(inf, recipient_cell);
    }
    return n_kernel_calls;
}

//Conditional entry method.
unsigned long long int Local_spread::localKeeling(Node* inf, Cell* recipient_cell)
{
    //Starts out by testing if we enter the cell
    size_t n_sus = recipient_cell->get_n_sus();
//    double p_over = inf->get_parent()->get_p_over(recipient_cell);
    double p_over = calc_infection_p_over_USDOS2(inf, recipient_cell);
    unsigned long long int n_kernel_calls = 1;
    //p_cell = p of 'entering' the recipient cell if all nodes of recipient were at
    //the closest point between the two cells.
    double p_cell = 1 - std::pow((1 - p_over), n_sus);
    if(common::unif_rand() <= p_cell) //Enter cell.
    {
        bool overestimated_infection_occured = false;
        double P = 1 - std::pow((1 - p_over), n_sus + 1);
        double p_over_inv = 1 / (1 - p_over);
        //Go through each susceptible farm and check if it would become infected
        //if it was at the closest distance. The probability has to be adjusted for
        //the condition that at least one of these 'hypothetical' infections must occur
        //(since thats what 'entering the cell' means). As soon as this criterion is
        //fulfilled there is no need to adjust probabilities anymore so the flag is set
        //to false.
        std::vector<Node*> sus_v = recipient_cell->get_nodes(0);
        for(Node* sus : sus_v)
        {
            if(overestimated_infection_occured or p_over == 1.0)
                P = 1.0;
            else
                P = 1 - ((1 - P) * p_over_inv);
            //Check for infection if overestimated probability.
            double r = common::unif_rand();
            if(r < (p_over / P))
            {
                //If the node would have been at the shortest distance it would
                //have become infected. Now see if the node becomes infected using
                //its actual position. Also, set flag to true since one overestimated
                //infection has occured.
                overestimated_infection_occured = true;
                double p = (this->*p_function)(inf, sus);
                n_kernel_calls++;
                if(r < (p / P))
                    sus->become_exposed();
            }
        }
    }
    return n_kernel_calls;
}

//Simple conditional subsample method based on one infectious node and one recipient cell. Unused in the main paper, roughly equivalent to CE in performance.
unsigned long long int Local_spread::localBinomial(Node* inf, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 1;
    size_t n_sus = recipient_cell->get_n_sus();
    if(n_sus == 0)
    {
        return 0;
    }
    double p_over = inf->get_parent()->get_p_over(recipient_cell);
    size_t n_over = 0;
    p_over = this->generate_bin_rv_lookup(n_sus, p_over, n_over);
    //Draw the number of nodes that would get infected if they were at the closest distance
    //between the cells.
//    size_t n_over = common::draw_binom(n_sus, p_over);
    if(n_over > 0)
    {
        std::vector<Node*> sus_nodes = recipient_cell->get_nodes(0);
        std::vector<Node*> random_sample;
        common::random_unique(sus_nodes, n_over, random_sample);
        for(Node* node : random_sample)
        {
            double p = (this->*p_function)(inf, node);
            n_kernel_calls++;
            if(common::unif_rand() <= p/p_over)
            {
                node->become_exposed();
            }
        }
    }
    return n_kernel_calls;
}

//Conditional subsample method used in main paper.
unsigned long long int Local_spread::localBinomial_CtoC(Cell* inf_cell, Cell* recipient_cell)
{
    unsigned long long int n_kernel_calls = 1;
    size_t n_sus = recipient_cell->get_n_sus();
    if(n_sus == 0)
    {
        return 0;
    }
    //Probability of the most infectious node in c1 infecting the most susc. node in c2 at the shortest poss. distance.
    double p_over = inf_cell->get_p_over(recipient_cell);
    //Probability of one farm becomes infected by any of the infectious nodes in c1 using p_over
    size_t n_inf = inf_cell->get_n_inf();
    double cumulative_p_over = 1.0 - std::pow(1.0-p_over, n_inf);
    //Draw the number of nodes that would get infected using the cumulative p_over.
    size_t n_over = 0;
    cumulative_p_over = this->generate_bin_rv_lookup(n_sus, cumulative_p_over, n_over);
//    size_t n_over = common::draw_binom(n_sus, cumulative_p_over);
    if(n_over > 0)
    {
        std::vector<Node*> sus_nodes = recipient_cell->get_nodes(0);
        std::vector<Node*> random_sample;
        common::random_unique(sus_nodes, n_over, random_sample);
        std::vector<Node*> inf_nodes = inf_cell->get_nodes(2);
        for(Node* sus_node : random_sample)
        {
            double cumulative_p = 1.0;
            for(Node* inf_node : inf_nodes)
            {
                cumulative_p *= 1.0 - (this->*p_function)(inf_node, sus_node);
                n_kernel_calls++;
            }
            cumulative_p = 1.0 - cumulative_p; //Cumulative probability of at least one node in inf_cell infecting this sus_node.
            if(common::unif_rand() <= cumulative_p / cumulative_p_over)
            {
                sus_node->become_exposed();
            }
        }
    }
    return n_kernel_calls;
}
