# this library produces the plots for this software package
import matplotlib.pyplot as plt
import numpy as np


def implied_time_scale(dictionary):
    """
    This function plots the implied time scale as a function of time

    Inputs: A dictionary containing the elements necessary for the plot
    Outputs: An ITS plot
    """
    time = dictionary["x_axis"]
    data = dictionary["data"]
    x_label = dictionary["x_label"]
    y_label = dictionary["y_label"]
    leg = dictionary["legend"]
    save = dictionary["save"]
    ub = dictionary["ITS_ub"]

    for k in range(data.shape[1]):
        label_temp = "ITS" + str(data.shape[1] - k)
        plt.plot(time, data[0, k, :], label=label_temp)

    plt.xlim(0.0, ub)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    # plt.legend(loc="center right")
    plt.savefig(save)
    plt.close()


def kernel_generator_plot(dictionary, close=True):
    """
    This function plots the memory kernel as a function of time.

    Inputs:
    1. mem - the memory kernel as a function of time
             dimension: (rd, rd, tsteps)

    Outputs:
    1. rd^2 plots for the memory kernel dynamics
    """
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # reading from the dictionary
    data = dictionary["data"]
    clr = dictionary["color"]
    mkr = dictionary["mkr"]
    lbl = dictionary["lbl"]
    title = dictionary["title"]
    num = dictionary["num"]

    # number of rows (columns)
    num_rows = data.shape[0]
    tsteps = data.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    # set length & width
    width = 12
    length = 12
    lw = 2.0

    # set tick parameters
    lb = min(timearray)
    ub = max(timearray) + 1
    sp = int(max(timearray) / 2)
    xticks = np.arange(lb, ub, sp)

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):
            # useful variables
            subplot_number = (num_rows * i) + j + 1

            # create subplots
            plt.subplot(num_rows, num_rows, subplot_number)
            plt.xlim(0, len(data[0, 0, :]))

            if (i == (num_rows - 1)):
                plt.xticks(xticks)
                plt.xlabel("t (step)", fontsize=15)

                if (j == 0):
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             label=lbl)
                else:
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             label=lbl)
            else:
                # omit the x-axis labels
                plt.gca().xaxis.set_major_locator(plt.NullLocator())

                if (j == 0):
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             label=lbl)

                else:
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             label=lbl)

            if ((i == (num_rows - 1)) and (j == (num_rows - 1))):
                plt.legend(loc='center right', frameon=False)

    plt.tight_layout(pad=num)
    plt.savefig(title)
    if close is True:
        plt.close()


def plot(dictionary, close=True):
    """
    Plots data that is given in the form of a dictionary

    Outputs:
    1. rd^2 plots of the dynamics as a function of time
    """
    # reading from the dictionary
    data = dictionary["data"]
    clr = dictionary["color"]
    mkr = dictionary["mkr"]
    lbl = dictionary["lbl"]
    title = dictionary["title"]

    # markevery (for convenience)
    # if clr is "forestgreen":
    #     me = 15
    # else:
    #     me = 15

    tsteps = data.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # number of rows (columns)
    num_rows = data.shape[0]
    tsteps = data.shape[-1]

    # set length & width
    width = 12
    length = 12
    lw = 3.0

    # set tick parameters
    lb = min(timearray)
    ub = max(timearray) + 1
    sp = int(max(timearray) / 2)
    xticks = np.arange(lb, ub, sp)
    yticks = np.arange(0.0, 1.1, 0.5)

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1

            # create subplots
            plt.subplot(num_rows, num_rows, subplot_number)
            # plt.xlim(0, len(data[0, 0, :]))
            plt.xlim(0, tsteps + 0.05)
            plt.ylim(0, 1.1)

            if (i == (num_rows - 1)):
                plt.xlabel("t (step)", fontsize=15)
                plt.xticks(xticks)

                if (j == 0):
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=8,
                             label=lbl)
                    plt.yticks(yticks)
                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=8,
                             label=lbl)

            else:
                if (j == 0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=8,
                             label=lbl)
                    plt.yticks(yticks)
                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=8,
                             label=lbl)

            if ((i == (num_rows - 1)) and (j == (num_rows - 1))):
                plt.legend(loc='upper right', frameon=False)

    plt.savefig(title)
    if close is True:
        plt.close()


def plot_RMSEs(dictionary, ub=1500.0, close=True):
    """
    Plots RMSE data from a dictionary.
    Inputs:
    1. dictionary
    2. ub
    3. close
    """
    ub = 1202
    width, length = 8, 5
    plt.figure(0, figsize=(width, length))
    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15) 
    plt.xticks(np.arange(0.0, int(ub) + 1, 200))

    # read from the dictionary
    x_axis = dictionary["x_axis"]
    data = dictionary["data"]
    clr = dictionary["color"]
    #mkr = dictionary["mkr"]
    mkr = "-"
    x_title = dictionary["xtitle"]
    title = dictionary["title"]
    
    # plt.plot(x_axis, data*0.0, color="salmon")
    plt.plot(x_axis, data,
             mkr,
             c=clr,
             markersize=4.0,
             linewidth=2.0)

    y_ub = 0.08
    # plt.ylim(0.0, y_ub)
    plt.xlim(-100.0, ub)
    plt.xlabel(x_title)
    plt.ylabel("RMSE")
    title = "./output/" + title
    plt.savefig(title)

    if close is True:
        plt.close()


def plots_for_andres(dictionary, close=True):
    """
    Plots the 1,1 TPM matrix dynamics through time.

    Outputs:
    1. rd^2 plots of the dynamics as a function of time
    """
    data = dictionary["data"]
    clr = dictionary["color"]
    lbl = dictionary["lbl"]

    tsteps = data.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # number of rows (columns)
    num_rows = data.shape[0]
    tsteps = data.shape[-1]

    # set length & width
    width = 15
    length = 15
    lw = 2.0

    # set tick parameters
    lb = min(timearray)
    ub = max(timearray)
    sp = int(max(timearray) / 2)
    xticks = np.arange(0, 21, 10)
    # yticks = np.arange(0.0, 1.1, 0.5)
    #plt.xlim(0.0, 21)
    #plt.ylim(0.9, 1)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.xlabel("time (step)", fontsize=15)
    plt.xticks(xticks)
    plt.plot(timearray, data[0, 0, :], "-", markersize=2, c=clr, label=lbl)
    X, Y = np.loadtxt('T11data.txt', delimiter=',', unpack=True)
    plt.plot(X, Y, "o", c="firebrick")

    plt.legend(loc='upper right', frameon=False)
    plt.savefig("./output/log_many_state_ITS.pdf")

    if close is True:
        plt.close()


def MS_ITS_plots(x, y, close=True):
    """
    Plots the 1,1 TPM matrix dynamics through time.

    Outputs:
    1. rd^2 plots of the dynamics as a function of time
    """

    plt.rc('xtick', labelsize=12)
    plt.rc('ytick', labelsize=12)

    # set length & width
    width = 15
    length = 15
    lw = 2.0

    # set tick parameters
    #lb = min(timearray)
    #ub = max(timearray)
    #sp = int(max(timearray) / 2)
    xticks = np.arange(0, 51, 25)
    # yticks = np.arange(0.0, 1.1, 0.5)
    # plt.xlim(0, len(data[0, 0, :]))
    # plt.ylim(-1.8, -1.25)
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.xlabel("time (step)", fontsize=12)
    plt.xticks(xticks)
    y[0] = np.nan
    y[1] = np.nan
    plt.plot(x, y, "-")

    # plt.legend(loc='upper right', frameon=False)
    plt.savefig("./output/many_state_ITS.pdf")

    if close is True:
        plt.close()


def plotter(dictionary, t_max, close=True):
    """
    Plots data that is given in the form of a dictionary

    Outputs:
    1. rd^2 plots of the dynamics as a function of time
    """

    # reading from the dictionary
    data = dictionary["data"]
    clr = dictionary["color"]
    mkr = dictionary["mkr"]
    lbl = dictionary["lbl"]
    title = dictionary["title"]
    lt = dictionary["lt"]

    tsteps = data.shape[-1]
    timearray = np.arange(0, t_max+1, lt)

    #plt.rc('xtick', labelsize=15)
    #plt.rc('ytick', labelsize=15)

    # number of rows (columns)
    num_rows = data.shape[0]
    tsteps = data.shape[-1]

    # set length & width
    width = 12
    length = 12
    lw = 3.0

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # set tick parameters
    lb = 0.0
    ub = t_max
    sp = int(t_max) / 2
    xticks = np.arange(lb, ub, sp)
    yticks = np.arange(0.0, 1.1, 0.5)

    # plot the matrix elements
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            plt.xlim(0, t_max)
            plt.ylim(0.0, 1.1)

            # create subplots
            plt.subplot(num_rows, num_rows, subplot_number)

            if (i == (num_rows - 1)):
                plt.xlabel("t (step)", fontsize=15)
                plt.xticks(xticks)

                if (j == 0):
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=4,
                             label=lbl)
                    plt.xlim(0, t_max)
                    plt.ylim(0.0, 1.1)
                    plt.yticks(yticks)

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=4,
                             label=lbl)
                    plt.xlim(0, t_max)
                    plt.ylim(0.0, 1.1)

            else:
                if (j == 0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=4,
                             label=lbl)
                    plt.yticks(yticks)
                    plt.xlim(0, t_max)
                    plt.ylim(0.0, 1.1)

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.plot(timearray, data[i, j, :],
                             mkr,
                             c=clr,
                             linewidth=lw,
                             markersize=4,
                             label=lbl)
                    plt.xlim(0, t_max)
                    plt.ylim(0.0, 1.1)

            if ((i == (num_rows - num_rows)) and (j == (num_rows - 1))):
                plt.legend(loc='lower right', frameon=False)

    plt.xlim(0, t_max)
    plt.ylim(0.0, 1.1)

    plt.savefig(title)
    #plt.show()
    if close is True:
        plt.close()



def plot_lit(dictionary, num, lt, lbl, show, t_max, close=False, save=False):
    """
    Plots data that is given in the form of a dictionary

    Outputs:
    1. rd^2 plots of the dynamics as a function of time
    """

    # reading from the dictionary
    data = dictionary["data"]
    clr = dictionary["color"]
    mkr = dictionary["mkr"]
    title = dictionary["title"]

    tsteps = data.shape[-1]
    timearray = np.arange(0, t_max + 1, lt)

    plt.rc('xtick', labelsize=15)
    plt.rc('ytick', labelsize=15)

    # number of rows (columns)
    num_rows = data.shape[0]
    tsteps = data.shape[-1]

    # set length & width
    width = 12
    length = 12
    lw = 3.0

    # define figure dimensions
    plt.figure(0, figsize=(width, length))
    plt.gca().yaxis.set_major_locator(plt.NullLocator())
    plt.gca().xaxis.set_major_locator(plt.NullLocator())
    plt.plot(timearray, data[num, num, :],
                             mkr,
                             linewidth=lw,
                             markersize=5,
                             label=lbl)
    
    if (show == 1):
        X, Y = np.loadtxt('./a2p_data/T44data.txt', delimiter=',', unpack=True)
        plt.plot(X, Y, "o", c="black")

    elif (show == 2):
        X, Y = np.loadtxt('./a2p_data/T22data.txt', delimiter=',', unpack=True)
        plt.plot(X, Y, "o", c="black")

    elif (show == 3):
        X, Y = np.loadtxt('./a2p_data/T33data.txt', delimiter=',', unpack=True)
        plt.plot(X, Y, "o", c="black")

    elif (show == 4):
        X, Y = np.loadtxt('./a2p_data/T11data.txt', delimiter=',', unpack=True)
        plt.plot(X, Y, "o", c="black")

    xticks = np.arange(0.0, t_max, 10.0)
    yticks = np.arange(0.6, 1.1, 0.5)
    plt.xlim(0.0, 100.5)
    plt.ylim(0.0, 1.1)
    plt.xticks(xticks)
    plt.yticks(yticks)
    plt.legend(loc='upper right', frameon=False)

    if save is True: 
        plt.savefig("./output/lit_comparison.pdf")

    # plt.show()
    if close is True:
        plt.close()
