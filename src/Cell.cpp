#include "Cell.h"
#include "common_functions.h"
#include "Config.h"
#include "Grid.h"
#include <algorithm>
#include <math.h>

unsigned int Cell::cells_created = 0; //Keeps track of total number of cells.

Cell::Cell(double lower_left_x, double lower_left_y, double width, double height,
           Grid* parent, const Config& C, bool divide_along_x) :
    width(width),
    height(height),
    //Corners
    ll(lower_left_x, lower_left_y),
    lr(lower_left_x + width, lower_left_y),
    ul(lower_left_x, lower_left_y + height),
    ur(lower_left_x + width, lower_left_y + height),
    //Grid object that the cell belongs to
    parent(parent),
    //Access to config options
    C(C),
    divide_along_x(divide_along_x)
{
    corners = {ll, lr, ul, ur};
    x_bounds.resize(2);
    y_bounds.resize(2);
    x_bounds[0] = lower_left_x;
    x_bounds[1] = lower_left_x + width;
    y_bounds[0] = lower_left_y;
    y_bounds[1] = lower_left_y + height;

    if(width == height)
    {
        quadratic = true;
    }
    id = cells_created;
    cells_created++;
}

Cell::Cell(double lower_left_x, double lower_left_y, double side,
           Grid* parent, const Config& C) :
    Cell(lower_left_x, lower_left_y, side, side, parent, C, false)
{

}

Cell::~Cell()
{
}

unsigned int Cell::get_n_cells()
{
    return cells_created;
}

unsigned int Cell::get_id()
{
    return id;
}

double Cell::get_max_susceptibility()
{
    return max_susceptibility;
}

double Cell::get_max_infectiousness()
{
    return max_infectiousness;
}

std::vector<Node*> Cell::get_nodes()
{
    return nodes;
}

//Get nodes of specific status.
std::vector<Node*> Cell::get_nodes(unsigned int status)
{
    std::vector<Node*>* vp;
    if(status == 0)
        vp = &sus_nodes;
    else if(status == 1)
        vp = &exp_nodes;
    else if(status == 2)
        vp = &inf_nodes;
    else
    {
        std::cout << "Unknown status in arg to Cell::get_nodes(arg). arg is " << status
                  << "Exiting" << std::endl;
        exit(EXIT_FAILURE);
    }
    return *vp;
}

std::vector<Node*> Cell::get_transmitting_nodes()
{
    return tra_nodes;
}

unsigned int Cell::get_n_nodes()
{
    return n_nodes;
}

size_t Cell::get_distance_index()
{
    return distance_index;
}

double Cell::get_shortestDistanceTo(Cell* c)
{
    return this->parent->get_cell_cell_distance(this->distance_index, c->get_distance_index());
}

double Cell::get_p_over(Cell* c)
{
    return this->parent->get_cell_cell_pover(this->distance_index, c->get_distance_index());
}

void Cell::init_node_statuses()
{
    sus_nodes.clear();
    exp_nodes.clear();
    inf_nodes.clear();
    det_nodes.clear();
    rem_nodes.clear();
    tra_nodes.clear();
    sus_nodes.reserve(n_nodes);
    exp_nodes.reserve(n_nodes);
    inf_nodes.reserve(n_nodes);
    det_nodes.reserve(n_nodes);
    rem_nodes.reserve(n_nodes);
    tra_nodes.reserve(n_nodes);

    //Go through each node and update its status and add it to the correct status vector.
    for(Node* n : nodes)
    {
        switch(n->get_status())
        {
        case(0):
            sus_nodes.push_back(n);
            break;
        case(1):
            exp_nodes.push_back(n);
            break;
        case(2):
            inf_nodes.push_back(n);
            break;
        case(3):
            det_nodes.push_back(n);
            break;
        case(4):
            rem_nodes.push_back(n);
        }
    }

    n_sus = sus_nodes.size();
    n_inf = inf_nodes.size() + det_nodes.size(); //Detected nodes are still infectious.
    //Fill a vector with the nodes that are transmitting infection (inf + det).
    tra_nodes.insert(tra_nodes.end(), inf_nodes.begin(), inf_nodes.end());
    tra_nodes.insert(tra_nodes.end(), det_nodes.begin(), det_nodes.end());
}

void Cell::init_transmission_parameters()
{
    //Saves the highest sus/trans of the nodes in this cell so that
    //they can be easily accessed when doing 'enter cell'-probabilities.
    for(Node* node : nodes)
    {
        if(node->get_susceptibility() > this->max_susceptibility)
        {
            this->max_susceptibility = node->get_susceptibility();
        }
        if(node->get_infectiousness() > this->max_infectiousness)
        {
            this->max_infectiousness = node->get_infectiousness();
        }
    }
}

//Resets the status of every node in the cell.
void Cell::reset()
{
    for(Node* n : nodes)
    {
        n->reset();
    }
    this->init_node_statuses(); //Reset all status vectors.
}

//After transmission has been evaluated, update each nodes' status and add them to the correct vector.
void Cell::commitTimeStep()
{
    sus_nodes.clear();
    exp_nodes.clear();
    inf_nodes.clear();
    det_nodes.clear();
    rem_nodes.clear();
    tra_nodes.clear();

    //Go through each node and update its status and add it to the correct status vector.
    for(Node* n : nodes)
    {
        n->update();
        switch(n->get_status())
        {
        case(0):
            sus_nodes.push_back(n);
            break;
        case(1):
            exp_nodes.push_back(n);
            break;
        case(2):
            inf_nodes.push_back(n);
            break;
        case(3):
            det_nodes.push_back(n);
            break;
        case(4):
            rem_nodes.push_back(n);
        }
    }

    n_sus = sus_nodes.size();
    n_inf = inf_nodes.size() + det_nodes.size(); //Detected nodes are still infectious.
    //Fill a vector with the nodes that are transmitting infection (inf + det).
    tra_nodes.insert(tra_nodes.end(), inf_nodes.begin(), inf_nodes.end());
    tra_nodes.insert(tra_nodes.end(), det_nodes.begin(), det_nodes.end());
}

void Cell::current_status(std::vector<size_t>& status_vec)
{
    status_vec[0] += sus_nodes.size();
    status_vec[1] += exp_nodes.size();
    status_vec[2] += inf_nodes.size();
    status_vec[3] += det_nodes.size();
    status_vec[4] += rem_nodes.size();
}

size_t Cell::get_n_sus()
{
    return n_sus;
}

size_t Cell::get_n_inf()
{
    return n_inf;
}

size_t Cell::get_time()
{
    return parent->get_time();
}

//Is a specific point on the xy plane inside this cell?
bool Cell::isWithin(Point p)
{
    if(p.get_x() >= x_bounds[0] and p.get_x() < x_bounds[1] and
       p.get_y() >= y_bounds[0] and p.get_y() < y_bounds[1])
    {
        return true;
    }
    else
    {
        return false;
    }
}

//Is a specific node within this cell.
bool Cell::isWithin(Node* n)
{
    return isWithin(n->get_position());
}

bool Cell::isWithinHRange(Point p)
{
    if(p.get_x() >= this->ll.get_x() && p.get_x() <= this->lr.get_x())
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool Cell::isWithinVRange(Point p)
{
    if(p.get_y() >= this->ll.get_y() && p.get_y() <= this->ul.get_y())
    {
        return true;
    }
    else
    {
        return false;
    }
}

std::vector<Point> Cell::get_corners()
{
    return corners;
}

void Cell::set_distance_index(size_t id)
{
    distance_index = id;
}

//Calculates the shortest distance to another cell when the grid is regular (static).
double Cell::calc_shortestDistanceSimple(Cell* other_cell)
{
    double shortest_distance = 999999999999999999999.0;
    std::vector<Point> other_corners = other_cell->get_corners();
    for(Point own_p : corners)
    {
        for(Point other_p : other_corners)
        {
            double current_distance = C.d_function(own_p, other_p);
            if(current_distance < shortest_distance)
            {
                shortest_distance = current_distance;
            }
        }
    }
    return shortest_distance;
}

//Calculates the shortest distance to another cell inside a grid with heterogeneously sized cells as created by the adaptive gridding method.
double Cell::calc_shortestDistanceComplex(Cell* other_cell)
{
    bool shortest_distance_found = false;
    std::vector<double> new_distances;
    new_distances.reserve(8);

    std::vector<Point> other_corners = other_cell->get_corners();

    for(Point other_p : other_corners)
    {
        for(Point own_p : corners)
        {
            //if c1 & c2 share a point they are adjacent and d = 0
            if(own_p == other_p)
            {
                new_distances.push_back(0.0);
                shortest_distance_found = true;
                break;
            }
        }
    }

    //If an edge of c2 overlaps with an edge in c1, the two points
    //between which the shortest distance is found is on a straight horizontal or
    //vertical line (in which case we don't need to use the distance function).
    //Check for horizontal alignment first.
    if(!shortest_distance_found)
    {
        for(Point other_p : other_corners)
        {
            if(this->isWithinHRange(other_p))
            {
                new_distances.push_back(std::abs(this->ll.get_y() - other_p.get_y()));
                new_distances.push_back(std::abs(this->ul.get_y() - other_p.get_y()));
                shortest_distance_found = true;
            }
        }
        for(Point own_p : this->corners)
        {
            if(other_cell->isWithinHRange(own_p))
            {
                new_distances.push_back(std::abs(other_corners[0].get_y() - own_p.get_y()));
                new_distances.push_back(std::abs(other_corners[2].get_y() - own_p.get_y()));
                shortest_distance_found = true;
            }
        }

        //Check for vertical alignment
        for(Point other_p : other_corners)
        {
            if(this->isWithinVRange(other_p))
            {
                new_distances.push_back(std::abs(this->ll.get_x() - other_p.get_x()));
                new_distances.push_back(std::abs(this->lr.get_x() - other_p.get_x()));
                shortest_distance_found = true;
            }
        }
        for(Point own_p : this->corners)
        {
            if(other_cell->isWithinVRange(own_p))
            {
                new_distances.push_back(std::abs(other_corners[0].get_x() - own_p.get_x()));
                new_distances.push_back(std::abs(other_corners[1].get_x() - own_p.get_x()));
                shortest_distance_found = true;
            }
        }
    }
    //No vertical or horizontal overlap, must use euclidean distance
    if(!shortest_distance_found)
    {
        for(Point own_p : this->corners)
        {
            for(Point other_p : other_corners)
            {
                new_distances.push_back(C.d_function(own_p, other_p));
            }
        }
    }

    //Get minimum from the vector and return it.
    return *std::min_element(new_distances.begin(), new_distances.end());
}

void Cell::addNode(Node* node)
{
    nodes.push_back(node);
    node->set_parent(this);
    n_nodes++;
}

void Cell::addNodes(std::vector<Node*> in_nodes)
{
    for(Node* n : in_nodes)
    {
        addNode(n);
    }
}

std::vector<Cell*> Cell::subdivideQuad(unsigned int max_nodes, bool force_division)
{
    //Tell a cell to subdivide into 4 smaller cells if number of nodes
    //is above a certain value. Only allowed for quadratic cells.
    //If force_division is set to true any conditions that control if division
    //will happen or not will be ignored and the cell will always be divided.
    if(!quadratic)
    {
        std::cout << "Attempt do subdivide a cell that is not quadratic into 4 children. "
                  << "Height, width: " << height << ", " << width << std::endl;
        exit(EXIT_FAILURE);
    }

    std::vector<Cell*> completed_cells;
    //If cell needs to subdivide, create 4 new cells with half the side
    if(nodes.size() > max_nodes or force_division)
    {
        std::vector<Cell*> children;
        double half_side = width * 0.5;
        children.push_back(new Cell(ll.get_x(), ll.get_y(), half_side, parent, C));
        children.push_back(new Cell(ll.get_x()+half_side, ll.get_y(), half_side, parent, C));
        children.push_back(new Cell(ll.get_x(), ll.get_y()+half_side, half_side, parent, C));
        children.push_back(new Cell(ll.get_x()+half_side, ll.get_y()+half_side, half_side, parent, C));


        std::vector<std::vector<Node*>> cell_nodes(4, std::vector<Node*>());
//      OLD  std::vector<double> log_c_sizes(4,0);
        std::vector<double> log_c_sizes;
        std::vector<double> c_sizes;
        for(size_t i = 0; i < children.size(); i++)
        {
            for(Node* node : nodes)
            {
                if(children[i]->isWithin(node))
                    cell_nodes[i].push_back(node);
            }
            //Count the number of nodes that will end up in each cell.
            //Only count variance for cells that have nodes
            if(cell_nodes[i].size() > 0)
                log_c_sizes.push_back(std::log(cell_nodes[i].size()));
                c_sizes.push_back(cell_nodes[i].size());
        }

        //Calculate new and original variance from the target value
        double c_var = common::variance_from(log_c_sizes, std::log(max_nodes));
        double p_var = double(std::log(this->get_n_nodes()) - std::log(max_nodes)) * double(std::log(this->get_n_nodes()) - std::log(max_nodes));

        if(c_var < p_var or force_division)
        {
            //Actually add the nodes to the cell.
            for(size_t i = 0; i < children.size(); i++)
            {
                for(Node* node : cell_nodes[i])
                {
                    children[i]->addNode(node);
                }
            }
            std::vector<Cell*> children_with_nodes;
            for(Cell* child : children)
            {
                if(child->get_n_nodes() > 0)
                {
                    children_with_nodes.push_back(child);
                }
                else
                {
                    delete child;
                }
            }
            //Subdivide the children that has been created at this level
            for(Cell* child : children_with_nodes)
            {
                std::vector<Cell*> grandchildren = child->subdivideQuad(max_nodes);
                //Add the new children to the complete list of cells.
                completed_cells.insert(completed_cells.end(),
                                       grandchildren.begin(), grandchildren.end());
            }
        }
        else
        {
            completed_cells.insert(completed_cells.end(),
                               this);
        }
    }
    else
    {
        completed_cells.insert(completed_cells.end(),
                               this);
    }
    return completed_cells;
}

std::string Cell::as_str()
{
    std::stringstream ss;
    ss << "Cell " << id << ": "
       << "ll (" << ll.get_x() << "; " << ll.get_y() << "); "
       << "lr (" << lr.get_x() << "; " << lr.get_y() << "); "
       << "ul (" << ul.get_x() << "; " << ul.get_y() << "); "
       << "ur (" << ur.get_x() << "; " << ur.get_y() << "); "
       << "N nodes: " << n_nodes;

    return ss.str();
}
