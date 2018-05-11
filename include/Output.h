#ifndef OUTPUT_H_INCLUDED
#define OUTPUT_H_INCLUDED

#include <string>
#include <fstream>
#include <vector>
#include <unordered_map>
#include <set>
#include "Grid.h"
#include "Node.h"

///Class to contain functions that handle the output from the simulations.

namespace output
{
    int do_header(std::string res_fn)
    {
        std::ofstream res_f(res_fn, std::ifstream::out);
        if(res_f.is_open())
        {
            res_f << "replicate\tseeded_node\tcum_n_inf\tmax_simul_n_inf\ttimestep_at_exit\t"
                     "n_kernel_calls\trep_time_ms\tinit_time_ms\tcell_size_mean(n_nodes)\t"
                     "cell_size_stdev_from_mean(n_nodes)\tcell_size_stdev_from_maxsize(nnodes)\t"
                     "n_cells\n";
            res_f.close();
            return 0;
        }
        else
        {
            std::cout << "Result file failed to open. Exiting..."
                      << std::endl; exit(EXIT_FAILURE);
        }
        return 1;
    }

    int do_group_header(std::string res_fn, std::vector<unsigned int> groups)
    {
        std::ofstream res_f(res_fn, std::fstream::out);
        if(res_f.is_open())
        {
            res_f << "Rep\tSeededGroups\t";
            for(int i=0; i<groups.size()-2; i++)
            {
                res_f << groups.at(i) << "\t";
            }
            res_f << groups.at(groups.size()-1) << std::endl;
            res_f.close();
            return 0;
        }
        else
        {
            std::cout << "Result file failed to open. Exiting..."
                      << std::endl; exit(EXIT_FAILURE);
        }
        return 1;
    }

    int do_output(std::string res_fn, Grid& G, Config& C, size_t time_step, unsigned int rep_i,
                          unsigned int rep_duration, unsigned int init_duration)
    {
        std::vector<size_t> current_status = G.get_current_status();
        std::vector<unsigned long long int> summary_vec = {current_status[1] + current_status[2] + current_status[3] + current_status[4],
                                                           G.get_max_n_inf(), time_step, G.get_replicate_kc()};

        std::ofstream res_f(res_fn, std::ifstream::app);
        if(res_f.is_open())
        {
            res_f << rep_i << "\t";
            res_f << "(";
            std::vector<Node*> seeded_nodes = G.get_seeded_nodes();
            for(size_t i=0; i<seeded_nodes.size()-1; i++)
            {
                res_f << seeded_nodes[i]->get_id() << ", ";
            }
            res_f << seeded_nodes.back()->get_id() << ")\t";

            for(size_t i=0; i<summary_vec.size(); i++)
            {
                res_f << summary_vec[i] << "\t";
            }
            res_f << rep_duration << "\t";
            res_f << init_duration << "\t";
            double cell_size_mean = G.get_cell_size_mean();
            res_f << cell_size_mean << "\t";
            res_f << G.get_cell_size_stdev_from(cell_size_mean) << "\t";
            res_f << G.get_cell_size_stdev_from(C.cell_max_farms) << "\t";
            res_f << G.get_n_cells() << "\n";
            res_f.close();
            return 0;
        }
        else
        {
            std::cout << "Result file not opened. Exiting..."
                      << std::endl; exit(EXIT_FAILURE);
        }
        return 1;
    }

    int do_group_output(std::string res_fn, unsigned int rep,
                        std::unordered_map<unsigned int, std::vector<Node*>> group_infections,
                        std::vector<unsigned int> groups, std::vector<Node*> seeded_nodes)
    {
        std::ofstream res_f(res_fn, std::fstream::app);
        if(res_f.is_open())
        {
            res_f << rep << "\t(";
            std::set<unsigned int> seeded_groups_set;
            for(Node* n : seeded_nodes)
            {
                seeded_groups_set.insert(n->get_group());
            }
            std::vector<unsigned int> seeded_groups;
            seeded_groups.assign(seeded_groups_set.begin(), seeded_groups_set.end());
            for(int i=0; i<seeded_groups.size()-1; i++)
            {
                res_f << seeded_groups.at(i) << "\t";
            }
            res_f << seeded_groups.back() << ")\t";
            for(int i=0; i<groups.size()-1; i++)
            {
                unsigned int group = groups.at(i);
                if(group_infections.find(group) == group_infections.end())
                {
                    res_f << 0 << "\t";
                }
                else
                {
                    res_f << group_infections.at(group).size() << "\t";
                }

            }
            unsigned int group = groups.at(groups.size()-1);
            if(group_infections.find(group) == group_infections.end())
            {
                res_f << 0 << std::endl;
            }
            else
            {
                res_f << group_infections.at(group).size() << std::endl;
            }

            res_f.close();
            return 0;
        }
        else
        {
            std::cout << "Result file failed to open. Exiting..."
                      << std::endl; exit(EXIT_FAILURE);
        }
        return 1;
    }

    int do_inf_coordinates_output(std::string inf_fn, Grid& G)
    {
        std::ofstream out_inf_coords(inf_fn);
        std::vector<Node*> nodes = G.get_all_nodes();
        for(Node* n : nodes)
        {
            out_inf_coords << n->get_position().get_x() << "\t"
                           << n->get_position().get_y() << "\t"
                           << n->get_time_of_infection() << std::endl;
        }
    }
}


#endif // OUTPUT_H_INCLUDED
