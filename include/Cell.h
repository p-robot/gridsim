#ifndef CELL_H
#define CELL_H

#include <vector>
#include <unordered_map>
#include "Node.h"
#include "Point.h"

class Grid;
class Config;
class Bin5tbl;

/**
Represents one cell inside the grid.
*/
class Cell
{
    public:
        /**
        Constructed using the lower left x & y coordinates as well
        as side length.
        */
        Cell();
        Cell(double lower_left_x, double lower_left_y, double side, Grid* parent, const Config& C);
        Cell(double lower_left_x, double lower_left_y, double sideH, double sideV, Grid* parent, const Config& C, bool divide_along_x);
        ~Cell();
        static unsigned int get_n_cells();
        unsigned int get_id(); ///Returns the cell's id.
        double get_max_susceptibility(); ///Returns the largest normalized susceptibility value of all member nodes.
        double get_max_infectiousness(); ///Returns the largest normalized infectiousness value of all member nodes.
        std::vector<Node*> get_nodes(); ///Returns all nodes.
        std::vector<Node*> get_nodes(unsigned int status); ///Returns nodes of specified status.
        std::vector<Node*> get_transmitting_nodes(); ///Returns all nodes that can transmit (status 2 & 3).
        unsigned int get_n_nodes(); ///Returns the number of nodes that belong to this cell.
        size_t get_distance_index(); ///Returns an index to where this cell is represented in distance vectors.
        double get_shortestDistanceTo(Cell* c); ///Returns the shortest distance between this cell and argument cell c.
        double get_p_over(Cell* c);
        void commitTimeStep(); ///Update all nodes statuses and put them in the correct status vectors.
        void current_status(std::vector<size_t>& status_vec); ///Adds n farms of different status to corresponding element in status_vec.
        size_t get_n_sus(); ///Returns the number of susceptible farms at the moment
        size_t get_n_inf(); ///Returns the number of infectious farms at the moment
        size_t get_time(); ///Gets the current timestep from the parent grid.
        bool isWithin(Point p); ///Checks if the point is within this cell.
        bool isWithin(Node* n); ///Checks if the node is within this cell.
        bool isWithinHRange(Point p); ///Checks if a point is within the horizontal range of the cell.
        bool isWithinVRange(Point p); ///Checks if a point is within the vertical range of the cell.
        std::vector<Point> get_corners(); ///Returns the corners in a vector (ll, lr, ul, ur)
        void set_distance_index(size_t id); ///Assigns id as an identifier to use when getting and setting distances.
        double calc_shortestDistanceSimple(Cell* other_cell); ///Calculate and return shortest distance from this cell to other for the static gridding method.
        double calc_shortestDistanceComplex(Cell* other_cell); ///Calculate and return shortest distance from this cell to other for the dynamic gridding method.
        void addNode(Node* node); ///Adds node to this cell and set the node's parent to this.
        void addNodes(std::vector<Node*> in_nodes); ///Adds nodes to this cell and set the nodes' parent to this.
        std::vector<Cell*> subdivideQuad(unsigned int max_nodes, bool force_division = false); ///Used to recursively create cells with dynamic size.
        void init_node_statuses(); ///Puts the nodes in the correct status vectors without updating them.
        void init_transmission_parameters(); ///Initiates maximum inf and max sus for this cell based on its nodes.
        void reset(); ///Incomplete
        std::string as_str(); ///Returns a summary of cell in string form.

    private:
        static unsigned int cells_created;
        Grid* parent; ///Pointer back up to the grid the cell exists within.
        unsigned int id;
        unsigned int n_nodes = 0;
        double width; ///Horizontal size of this cell.
        double height; ///Vertical size of cell.
        size_t n_sus = 0;
        size_t n_inf = 0;
        Point ll, lr, ul, ur; ///The four corners, lower left, lower right, upper left, upper right.
        const Config& C;
        bool divide_along_x = true;
        double max_susceptibility = 0.0;
        double max_infectiousness = 0.0;
        size_t distance_index = 0;
        bool quadratic = false;
        std::vector<Point> corners; ///Corners in a vector (lower left, lower right, upper left, upper right)
        std::vector<double> x_bounds, y_bounds; ///Max/min coordinates of this cells x/y dimensions.
        std::vector<Node*> nodes; ///All nodes that belong to this cell.
        std::vector<Node*> sus_nodes; ///All susceptible nodes in this cell.
        std::vector<Node*> exp_nodes; ///All exposed nodes in this cell.
        std::vector<Node*> inf_nodes; ///All infectious nodes in this cell.
        std::vector<Node*> det_nodes; ///All detected nodes in this cell.
        std::vector<Node*> tra_nodes; ///All currently transmitting nodes.
        std::vector<Node*> rem_nodes; ///All removed nodes in this cell.

        k_fun_ptr k_function; ///Pointer to kernel function in use.
        d_fun_ptr d_function; ///Pointer to distance function in use.
};

#endif // CELL_H
