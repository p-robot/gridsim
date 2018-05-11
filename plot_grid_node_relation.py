import numpy as np
from matplotlib import pyplot as plt
from matplotlib import colors
import matplotlib

def prettify_axes(ax):
    linewidth = 0.5

    #ax.tick_params(top='off', bottom = 'off', left = 'off', right='off', direction = 'out')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    
    ax.set_frame_on(False)
    ax.axis('equal')

def main():
    fontdict = {'size': 9}
    matplotlib.rc('font', **fontdict)

    fig, ax = plt.subplots(1,1)
    cmap = cmap = plt.get_cmap('winter')
    
    #Grid
    cells = list()
    cell_colors = list()
    fg = open("_grid.dat", 'r')
    for line in fg.readlines():
        l = line.strip()
        l = l.split('\t')
        coords = np.zeros((2,5))
        i = 0
        for x, y in zip(l[0::2], l[1::2]):
            coords[0,i] = x
            coords[1,i] = y
            i += 1
        coords[0,4] = coords[0,0]
        coords[1,4] = coords[1,0]
        #Switch place between UL and UR
        tempcULx = coords[0,2]
        tempcULy = coords[1,2]
        coords[0,2] = coords[0,3]
        coords[1,2] = coords[1,3]
        coords[0,3] = tempcULx
        coords[1,3] = tempcULy
        cells.append(coords)
        cell_colors.append(int(l[8]))
    
    norm = colors.Normalize(min(cell_colors), max(cell_colors))
    
    for i in range(len(cells)):
        ax.plot(cells[i][0,:], cells[i][1,:], lw = 0.25, c = cmap(norm(cell_colors[i])))
        # ax.plot(cells[i][0,:], cells[i][1,:], lw = 1.5, c = 'k')
    fg.close()
    
    #Nodes
    fn = open("_nodes.dat", 'r')
    n_coords = [[],[]]
    color_list = []
    
    for line in fn.readlines():
        l = line.strip()
        l = l.split('\t')
        n_coords[0].append(float(l[0]))
        n_coords[1].append(float(l[1]))
        color_list.append(float(l[2]))

    x_a = np.array(n_coords[0])
    y_a = np.array(n_coords[1])
    c_a = np.array(color_list)

    ax.scatter(x_a, y_a, c = c_a, norm = norm, cmap = cmap, s = .5, linewidth = 0)
    prettify_axes(ax)
    # plt.show()
    fig.set_size_inches(11.7*0.9, 8.27*0.9)
    plt.tight_layout()
    plt.savefig('node-grid-relationship.png', dpi = 1200)

if __name__ == '__main__':
    main()