#include "Shipment_handler.h"
#include "Node.h"
#include "common_functions.h"

Shipment_handler::Shipment_handler(std::string shipment_fn,
                                   std::unordered_map<unsigned int, Node*>& node_map)
{
    std::cout << "Loading shipments..." << std::flush;
    TimePoint process_start = common::get_wall_time();

    std::ifstream shipment_file(shipment_fn);
    if(shipment_file.is_open())
    {
        common::skipBOM(shipment_file);
        unsigned int n_lines = common::get_n_lines(shipment_file);
        shipment_vector.reserve(n_lines);
        std::string line;
        unsigned int line_number = 0;
        while(std::getline(shipment_file, line))
        {
            std::vector<std::string> line_vector = common::split(line, '\t');
            if(!line_vector.empty() and line_vector.size() == 5)
            {
                Shipment* s = new Shipment;
                s->id = common::stringToNum<unsigned int>(line_vector[0]);
                unsigned int o_node_id = common::stringToNum<unsigned int>(line_vector[1]);
                s->origin = node_map.at(o_node_id);
                unsigned int d_node_id = common::stringToNum<unsigned int>(line_vector[2]);
                s->destination = node_map.at(d_node_id);
                s->timestep = common::stringToNum<unsigned int>(line_vector[3]);
                s->n = common::stringToNum<unsigned int>(line_vector[4]);
                shipment_vector.push_back(s);
                time_shipment_map[s->timestep].push_back(s);
            }
            else
            {
                std::cout << "Error found in shipment file at line " << line_number << ". (" << line << ")." << std::endl;
                std::cout << "Line should read shipment_id\tfrom_id\tto_id\timestep\num_animals (tab separated)." << std::endl;
                exit(EXIT_FAILURE);
            }
        }
    }
    else
    {
        std::cout << "Failed to open shipment file " << shipment_fn << "." << std::endl;
        exit(EXIT_FAILURE);
    }
    unsigned long long int process_duration = common::delta_duration(process_start);
    std::cout << "complete (" << process_duration << " ms)." << std::endl;
}

Shipment_handler::~Shipment_handler()
{
    for(Shipment* s : shipment_vector)
        delete s;
}

void Shipment_handler::make_transmission(unsigned int t)
{
    auto it = time_shipment_map.find(t);
    if(it != time_shipment_map.end())
    {
        std::vector<Shipment*> ship_v = (it->second);
        for(Shipment* s : ship_v)
        {
            unsigned int origin_status = s->origin->get_status();
            if(origin_status == 2 or origin_status == 3)
            {
                if(s->destination->get_status() == 1)
                {
                    s->destination->become_exposed();
                }
            }
        }
    }
}
