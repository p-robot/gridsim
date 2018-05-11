This is gridsim, a small kernel based spatially explicit livestock disease
simulation. Gridsim is primarily intended to be a reference to the methods
described in the main paper, and as an example of how to implement them. As
such, it was never intended to be a disease simulation tool and the authors
make no claims about its power as a predictive model. 

After compiling and assuming that the binary is called gridsim, run the
program as "gridsim <config_file>" from the command line. Provided together
with the code is an example config file which can be easily be modified to
suit your needs. 
If you want to use your own landscape file simply construct a text file with
the columns [farm id, group, x-coordinate, y-coordinate, number of animal 
species 1, number of animal species 2]. Each line is a node and the file should
not contain a header. Group must be a number that indicates that the farm
belongs to for instance a region such as a county or similar, although it is
not necessary and every node can belong to the same group.
