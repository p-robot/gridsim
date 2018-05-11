#ifndef SHIPMENT_HANDLER_H
#define SHIPMENT_HANDLER_H

#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <sstream>

class Node;
struct Shipment;

class Shipment_handler
{
    public:
        Shipment_handler(std::string shipment_fn, std::unordered_map<unsigned int, Node*>& node_map);
        ~Shipment_handler();
        void make_transmission(unsigned int t); ///Perform all shipment-related transmission for the given time step.
    private:
        std::vector<Shipment*> shipment_vector;
        std::unordered_map<unsigned int, std::vector<Shipment*>> time_shipment_map;
};

struct Shipment
{
    Node* origin;
    Node* destination;
    unsigned int id;
    unsigned int timestep;
    unsigned int n; //Number of animals.
};

#endif // SHIPMENT_HANDLER_H
