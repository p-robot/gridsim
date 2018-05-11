#include <iostream>
#include <vector>
#include <ctime>
#include <fstream>
#include <sstream>
#include "Config.h"
#include "Grid.h"
#include "Node.h"
#include "common_functions.h"
#include "Output.h"

int main(int argc, char* argv[])
{
    //Read config filenam as command line argument
    if(argc != 2)
    {
        std::cout << "usage: " << argv[0] << " <config>" << std::endl;
        exit(EXIT_FAILURE);
    }
    //Construc t config object based on config file
    Config C = Config(argv[1]);

    //Write headers to output files.
    for(int level : C.output_levels)
    {
        std::stringstream ss;
        ss << level;
        std::string res_fn = C.batch_name + "_" + ss.str() + "_res.txt";
        output::do_header(res_fn);
    }


    TimePoint init_start = common::get_wall_time();
    //Construct the grid
    Grid G = Grid(C);

    if(C.make_plots)
    {
        G.plotGridAndNodes(false);
        G.plotCellSizeHist(C.cell_max_farms);
    }

    unsigned long long int init_duration = common::delta_duration(init_start);

    //Results output at the time step of outbreak dying.
    std::string end_fn = C.batch_name + "_end_res.txt";
    output::do_header(end_fn);

    //Number of infected per replicate and group as defined by column 2 in the landscape file.
    std::string group_out_fn = C.batch_name + "_group.txt";
    std::vector<unsigned int> groups = G.get_groups();
    output::do_group_header(group_out_fn, groups);

    int n_replicates = C.n_replicates;
    //If seeed method is 4 one node in each group will be randomly seeded and n replicates will be = n groups.
    if(C.seed_method == 4)
    {
        n_replicates = G.get_n_groups();
    }

    //Each rep starts here
    for(unsigned int rep_i = 0; rep_i < n_replicates; rep_i++)
    {
        unsigned int time_step = 0;
        unsigned int stopping_criterion_index = 0;
        //Reset the all the nodes' statuses each rep.
        G.reset();
        //Seed according to seed method
        std::vector<Node*> seeded_nodes = G.seed();

        //Stringsream to store number in each class for every timestep.
        std::stringstream incidence_stream;
        incidence_stream << "t\tn_sus\tn_exp\tn_inf\tn_det\tn_rem" << std::endl;
        //Get current numbers in each class.
        std::vector<size_t> current_status = G.get_current_status();
        //Insert numbers into stream.
        incidence_stream << time_step << "\t" << current_status[0] << "\t" //n sus
                         << current_status[1] << "\t" //n exposed
                         << current_status[2] << "\t" //n infectious
                         << current_status[3] << "\t" //n detected
                         << current_status[4] << std::endl; //n removed

        std::cout << "Replicate " << rep_i << " starting. " << std::endl
                  << " Status is: susceptible " << current_status[0]
                  << "; exposed " << current_status[1] << "; infectious "
                  << current_status[2] << "; detected " << current_status[3]
                  << "; removed " << current_status[4] << std::endl;

        //Start disease simulation.
        int exit_status = 0;
        TimePoint rep_start = common::get_wall_time();
        time_step += 1;
        while(time_step <= C.max_timesteps and exit_status == 0)
        {
            //Move forward one timestep and perform all transmissions.
            G.incrementTime();
            unsigned long long int rep_duration = common::delta_duration(rep_start);

            //Get numbers in each model compartment after transmissions.
            current_status = G.get_current_status();
            unsigned int current_n_inf = current_status[1] + current_status[2] + current_status[3];
            unsigned int cum_infected_nodes = current_n_inf + current_status[4];

            //Add numbers in infection classes to stream.
            incidence_stream << time_step << "\t" << current_status[0] << "\t" //n sus
                             << current_status[1] << "\t" //n exposed
                             << current_status[2] << "\t" //n infectious
                             << current_status[3] << "\t" //n detected
                             << current_status[4] << std::endl; //n removed

            //Outbreak died out.
            if(current_n_inf  == 0)
            {
                exit_status = 1;
                break;
            }

            //Print results for the outbreak stages that have been reached.
            while(cum_infected_nodes >= C.output_levels.at(stopping_criterion_index))
            {
                std::stringstream ss;
                ss << C.output_levels.at(stopping_criterion_index);
                std::string res_fn = C.batch_name + "_" + ss.str() + "_res.txt";
                output::do_output(res_fn, G, C, time_step, rep_i, rep_duration, init_duration);
                if(stopping_criterion_index >= C.output_levels.size()-1)
                {
                    exit_status = 2;
                    break;
                }
                stopping_criterion_index++;
            }
            //Increse time step.
            time_step++;
        }

        //Rep finished.
        unsigned long long int rep_duration = common::delta_duration(rep_start);
        output::do_output(end_fn, G, C, time_step, rep_i, rep_duration, init_duration);
        std::string inf_fn = C.batch_name + "_inf_" + std::to_string(rep_i) + ".txt";
//        output::do_inf_coordinates_output(inf_fn, G);

        std::unordered_map<unsigned int, std::vector<Node*>> group_infections = G.get_inf_nodes_by_group();
        output::do_group_output(group_out_fn, rep_i, group_infections, groups, seeded_nodes);

        //Output the stream of number of cases in each class to file.
//        std::string incidence_fn = C.batch_name + "_incidence_" + std::to_string(rep_i) + ".txt";
//        std::ofstream incidence_fs(incidence_fn);
//        if(incidence_fs.is_open())
//        {
//            incidence_fs << incidence_stream.rdbuf();
//            incidence_fs.close();
//        }
//        else
//        {
//            std::cout << "Failed to open incidence file (" << incidence_fn << "). No incidence output made." << std::endl;
//        }


        if(exit_status == 1)
        {
            std::cout << "No more transmission sources. Early exit at timestep " << time_step
                      << " (replicate " << rep_i << ")." << std::endl;
        }
        else if(exit_status == 2)
        {
            std::cout << "Reached " << C.output_levels[C.output_levels.size()-1]
                      << " total infected nodes. Early exit." << std::endl;
        }
        else
        {
            std::cout << time_step-1 << " time steps reached. Exiting (replicate " << rep_i << ")." << std::endl;
        }

    }
    return 0;
}
