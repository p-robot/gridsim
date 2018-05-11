#ifndef NODE_H
#define NODE_H

#include "Point.h"
#include "Config.h"
class Cell;
class Grid;

/**
Describes one node in the landscape. Status is described
as an integer, 0 = susceptible, 1 = exposed, 2 = infectious,
3 = detected, 4 = removed.
*/
class Node
{
    public:
        /**
        Construct using x and y coordinates, number of animals, id and config object.
        */
        Node(double x, double y, unsigned int num_sp1, unsigned int num_sp2, unsigned int id,
             unsigned int group, const Config& C, Grid* G);
        ~Node();
        Point get_position(); ///Returns the position of this node. Inlined
        unsigned int get_num_sp1(); ///Returns the number of animals of species 1 of this node.
        unsigned int get_num_sp2(); ///Returns the number of animals of species 2 of this node.
        unsigned int get_id(); ///Returns the node's id.
        unsigned int get_group(); ///Returns the node's group.
        double get_susceptibility(); ///Returns node's normalized susceptibility.
        double get_infectiousness(); ///Returns node's normalized infectiousness.
        int get_status(); ///Returns the current status of this node.
        int get_time_of_infection(); ///Returns the time step that this node got infected. Returns -1 if not infected.
        Cell* get_parent(); ///Returns the parent cell of this node.
        void set_susceptibility(double s); ///Sets the total normalized susceptibility of this node.
        void set_infectiousness(double i); ///Sets the total normalized infectiousness of this node.
        void set_status(int val); ///Sets the infections status of this node.
        void set_parent(Cell* in_parent); ///Sets the parent cell of this noe to arg.
        void update(); ///Update status for this timestep.
        void reset(); ///Reset this node to its original state (susceptible).
        void become_seeded(); ///Seed this node.
        void become_exposed(); ///Expose this node to infection.
        void become_infectious(); ///Turns this node infectious.
        void become_detected(); ///Turns this node into detected.
        void become_removed(); ///Turns this node into removed.

    private:
        Point position; ///Position.
        unsigned int num_sp1 = 0; ///Number of animals of species 1 of this node.
        unsigned int num_sp2 = 0; ///Number of animals of species 2 of this node.
        unsigned int id; ///Identifier for this node.
        unsigned int group; ///Which group does this node belong to. Eg. county.
        const Config& C; ///Config parameters.
        Grid* parent_grid; ///Pointer back to the grid object that the node belongs to.
        double susceptibility = 1.0; ///Total normalized susceptibility of this node.
        double infectiousness = 1.0; ///Total normalized infectiousness of this node.
        double current_cell_p_over = 1.0;
        unsigned int status = 0; ///Current status of node.
        int infection_time = -1;
        int countdown = 0; ///Keeps track of when its time for transition.
        Cell* parent; ///Parent cell.
};

inline Point Node::get_position()
{
    return position;
}

#endif // NODE_H
