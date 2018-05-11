import numpy as np
from matplotlib import pyplot as plt
import matplotlib

def main():
    fig, ax = plt.subplots(1,1)
    
    sizes = list()
    f = open('_cellsize.dat', 'r')
    max_nodes = f.readline()
    for line in f.readlines():
        sizes.append(int(line.strip()))
    f.close()
    
    ax.hist(sizes)
    ax.axvline(float(max_nodes))
    plt.savefig('cell-size-histogram.png', dpi = 300)
    
if __name__ == '__main__':
    main()
    