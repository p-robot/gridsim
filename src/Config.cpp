#include <fstream>
#include <sstream>
#include <algorithm>

#include "common_functions.h"
#include "Config.h"

Config::Config(std::string config_file)
{
	// Read config file and store each paramter in parameter_vector.
	std::ifstream f(config_file);
	if(f.is_open()){
		while(!f.eof()){ // until end of file
			std::string line;
			std::stringstream ss;
			// First get an entire line from the config file
			std::getline(f, line);
			// Put that in a stringstream (required for getline) and use getline again
			// to only read up until the comment character (delimiter).
			ss << line;
			std::getline(ss, line, '#');
			// If there is whitespace between the value of interest and the comment character
			// this will be read as part of the value. Therefore strip it of whitespace first.
			line.erase(std::remove_if(line.begin(), line.end(), isspace), line.end());
			if(line.size() != 0) //If there was anything but whitespace to the left of the delimiter.
            {
                parameter_vector.emplace_back(line);
            }
		}
	}

    if(parameter_vector.size() != 22)
    {
        std::cout << "Error in config file. Should be 22 parameters, found "
                  << parameter_vector.size() << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << parameter_vector.size() << " parameters successfully read from configuration file." << std::endl;

    if(parameter_vector[0] != "*" and parameter_vector[0] != "")
        batch_name = parameter_vector[0];

    if(parameter_vector[1] != "*" and parameter_vector[1] != "")
        node_file = parameter_vector[1];

    if(parameter_vector[2] != "*" and parameter_vector[2] != "")
        max_timesteps = common::stringToNum<unsigned int>(parameter_vector[2]);

    if(parameter_vector[3] != "*" and parameter_vector[3] != "")
        n_replicates = common::stringToNum<unsigned int>(parameter_vector[3]);

    if(parameter_vector[4] != "*" and parameter_vector[4] != "")
        opt_method = common::stringToNum<unsigned int>(parameter_vector[4]);

    if(parameter_vector[5] != "*" and parameter_vector[5] != "")
        local_kernel = common::stringToNum<unsigned int>(parameter_vector[5]);

    if(parameter_vector[6] != "*" and parameter_vector[6] != "")
        lookup_table = common::stringToNum<unsigned int>(parameter_vector[6]);

    if(parameter_vector[7] != "*" and parameter_vector[7] != "")
        cell_method = common::stringToNum<unsigned int>(parameter_vector[7]);

    if(parameter_vector[8] != "*" and parameter_vector[8] != "")
        n_cells_1d = common::stringToNum<unsigned int>(parameter_vector[8]);

    if(parameter_vector[9] != "*" and parameter_vector[9] != "")
        cell_max_farms = common::stringToNum<unsigned int>(parameter_vector[9]);

    if(parameter_vector[10] != "*" and parameter_vector[10] != "")
        inc_time = common::stringToNum<unsigned int>(parameter_vector[10]);

    if(parameter_vector[11] != "*" and parameter_vector[11] != "")
        det_time = common::stringToNum<unsigned int>(parameter_vector[11]);

    if(parameter_vector[12] != "*" and parameter_vector[12] != "")
        rem_time = common::stringToNum<unsigned int>(parameter_vector[12]);

    if(parameter_vector[13] != "*" and parameter_vector[13] != "")
        seed_method = common::stringToNum<unsigned int>(parameter_vector[13]);

    if(parameter_vector[14] != "*" and parameter_vector[14] != "")
        seed_n = common::stringToNum<unsigned int>(parameter_vector[14]);

    if(parameter_vector[15] != "*" and parameter_vector[15] != "")
        seed_file = parameter_vector[15];

    if(parameter_vector[16] != "*" and parameter_vector[16] != "")
        output_levels = common::split_stoi(parameter_vector[16], ',');

    if(parameter_vector[17] != "*" and parameter_vector[17] != "")
    {
        if(parameter_vector[17].compare("0") == 0)
            make_plots = false;
        else
            make_plots = true;
    }

    if(parameter_vector[18] != "*" and parameter_vector[18] != "")
    {
        std::vector<double> sus_scale_v = common::split_stod(parameter_vector[18], ',');
        S1 = sus_scale_v[0];
        S2 = sus_scale_v[1];
    }

    if(parameter_vector[19] != "*" and parameter_vector[19] != "")
    {
        std::vector<double> inf_scale_v = common::split_stod(parameter_vector[19], ',');
        T1 = inf_scale_v[0];
        T2 = inf_scale_v[1];
    }

    if(parameter_vector[20] != "*" and parameter_vector[20] != "")
    {
        std::vector<double> sus_exp_v = common::split_stod(parameter_vector[20], ',');
        q1 = sus_exp_v[0];
        q2 = sus_exp_v[1];
    }

    if(parameter_vector[21] != "*" and parameter_vector[21] != "")
    {
        std::vector<double> inf_exp_v = common::split_stod(parameter_vector[21], ',');
        p1 = inf_exp_v[0];
        p2 = inf_exp_v[1];
    }

//    if(parameter_vector[22] != "*" and parameter_vector[21] != "")
//    {
//        shipment_file = parameter_vector[22];
//        shipments_on = true;
//    }

    //Assign pointers to functions
    if(local_kernel == 1)
        k_function = &common::USDOS2_kernel;
    else if(local_kernel == 2)
        k_function = &common::gaussian_kernel;
    else if(local_kernel == 3)
        k_function = &common::USDOS2_kernel_half_gamma;
    else if(local_kernel == 4)
        k_function = &common::USDOS2_kernel_double_gamma;
    else if(local_kernel == 5)
        k_function = &common::brand_kernel_a3;
    else if(local_kernel == 6)
        k_function = &common::brand_kernel_a4;
    else if(local_kernel == 7)
        k_function = &common::brand_kernel_a5;
    else if(local_kernel == 8)
        k_function = &common::hayama_kernel;
    else if(local_kernel == 9)
        k_function = &common::USDOS2_kernel_v2;

    if(distance_method == 2)
        d_function = &common::sq_distance;
    else
        d_function = &common::eucl_distance;

    //Set up other parameters.
    if(opt_method == 4)
    {
        opt_on_node_level = false;
    }
    else
    {
        opt_on_node_level = true;
    }

    if (seed_method == 2)
    {
        std::ifstream f(seed_file, std::ifstream::in);
        if(f.is_open())
        {
            seed_n = common::get_n_lines(f);
        }
        else
        {
            std::cout << "Failed to open seed file " << seed_file << "." << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    else if(seed_method == 3)
    {
        specific_seed_id = seed_n;
        seed_n = 1;
    }
    else if(seed_method == 4)
    {
        seed_n = 1;
    }
}

Config::~Config()
{
    //dtor
}

std::string Config::as_str()
{
    std::stringstream ss;
    ss << ";batch_name=" << batch_name
                        << ";node_file=" << node_file
                        << "max_timesteps=" << max_timesteps
                        << ";n_replicates=" << n_replicates
                        << ";opt_method=" << opt_method
                        << ";cell_method=" << cell_method
                        << ";n_cells_1d=" << n_cells_1d
                        << ";cell_max_farms=" << cell_max_farms
                        << ";inc_time=" << inc_time
                        << ";det_time=" << det_time
                        << ";rem_time=" << rem_time
                        << ";seed_method=" << seed_method
                        << ";seed_file=" << seed_file
                        << "\n";
    return ss.str();
}
