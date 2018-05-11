#include "Node.h"
#include "Cell.h"
#include "Grid.h"

Node::Node(double x, double y, unsigned int num_sp1,
           unsigned int num_sp2, unsigned int id,
           unsigned int group, const Config& C, Grid* G) :
    position(x, y), num_sp1(num_sp1), num_sp2(num_sp2),
    id(id), group(group), C(C), parent_grid(G)
{

}

Node::~Node()
{
    //dtor
}

unsigned int Node::get_num_sp1()
{
    return num_sp1;
}

unsigned int Node::get_num_sp2()
{
    return num_sp2;
}

unsigned int Node::get_id()
{
    return id;
}

unsigned int Node::get_group()
{
    return group;
}

double Node::get_susceptibility()
{
    return susceptibility;
}

double Node::get_infectiousness()
{
    return infectiousness;
}

int Node::get_status()
{
    return status;
}

int Node::get_time_of_infection()
{
    return infection_time;
}

Cell* Node::get_parent()
{
    return parent;
}

void Node::set_susceptibility(double s)
{
    susceptibility = s;
}

void Node::set_infectiousness(double i)
{
    infectiousness = i;
}

void Node::set_status(int val)
{
    status = val;
}

void Node::set_parent(Cell* in_parent)
{
    parent = in_parent;
}

//Move the infection process forward one timestep for this node.
void Node::update()
{
    switch(status)
    {
    case(1): //If current status is exposed
        countdown--;
        if(countdown < 0)
            this->become_infectious();
        break;
    case(2): //If current status is infectious
        countdown--;
        if(countdown < 1)
            this->become_removed();
        break;
    case(3): //If status is detected
        countdown--;
        if(countdown < 1)
            this->become_removed();
        break;
    default:
        break;
    }
}

void Node::reset()
{
    status = 0;
    countdown = 0;
    infection_time = -1;
}

void Node::become_seeded()
{
    become_exposed();
    this->parent->init_node_statuses();
}

void Node::become_exposed()
{
    status = 1;
    countdown = C.inc_time;
    infection_time = parent->get_time();
    parent_grid->add_to_group_inf(this, group);
}

void Node::become_infectious()
{
    status = 2;
    countdown = C.det_time;
}

void Node::become_detected()
{
    status = 3;
    countdown = C.rem_time;
}

void Node::become_removed()
{
    status = 4;
    countdown = 0;
}
