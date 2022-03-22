import matplotlib.pyplot as plt

def Plotter(fig, axs, xy) -> NULL:
    """Takes a Figure from matplotlib, the array for the figure, and a list of data to plot,
    data in xy stored as a list of lists where xy = [["title", "namex", "namey", [x], [y]], [...]] 
    where x and y can be lists of plots themselves
    """
    index = 0
    for i in axs.reshape(-1):
        i.set_title(xy[index][0])
        i.set_xlabel(xy[index][1])
        i.set_ylabel(xy[index][2])
        if isinstance(xy[index][3], list) and isinstance(xy[index][4], list):
            for (j, k) in zip(xy[index][3], xy[index][4]):
                i.plot(j, k)

        elif isinstance(xy[index][3], list):
            for j in xy[index][3]:
                i.plot(j, xy[index][4])

        elif isinstance(xy[index][4], list):
            for j in xy[index][4]:
                i.plot(xy[index][3], j)

        else:
            i.plot(xy[index][3], xy[index][4])


        index += 1
    plt.tight_layout()