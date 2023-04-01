import matplotlib.pyplot as plt
import numpy as np


def exact_plot(data1, marker="-", label1=""):
    """
    This function plots the matrix elements of the various forms of the dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 

    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0,len(data1[:,:,0]))
            plt.ylim(0,1.1)

            if (i == 2):
                plt.xlabel("t (step)",fontsize=15)
                plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray) / 2)))
                if (j == 0):
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray) / 2)))
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)

            else:
                if (j==0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else: 
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=3.0, markersize=5, markevery=500, label=label1)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)
    plt.show()


def plot_data(data1, data3, data4, data5, msm_axis, color="aqua", color2="red", color3="red", 
              marker="--", marker2="--", marker3="--", label1="", label3="", label4="", label5=""):
    """
    This function plots the matrix elements of the various forms of the dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 

    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0,len(data1[:,:,0]))
            plt.ylim(0,1.1)

            if (i == 2):
                plt.xlabel("t (step)",fontsize=15)
                plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray)) / 2))
                if (j == 0):
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data3[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label3)
                    plt.plot(timearray, data4[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label4)
                    plt.plot(msm_axis, data5[i, j, :], marker3, c=color3, linewidth=3.0, markersize=7, label=label5)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, 5000))
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data3[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label3)
                    plt.plot(timearray, data4[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label4)
                    plt.plot(msm_axis, data5[i, j, :], marker3, c=color3, linewidth=3.0, markersize=7, label=label5)

            else:
                if (j==0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data3[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label3)
                    plt.plot(timearray, data4[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label4)
                    plt.plot(msm_axis, data5[i, j, :], marker3, c=color3, linewidth=3.0, markersize=7, label=label5)
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else: 
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=3.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data3[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label3)
                    plt.plot(timearray, data4[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label4)
                    plt.plot(msm_axis, data5[i, j, :], marker3, c=color3, linewidth=3.0, markersize=7, label=label5)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)

    plt.savefig('dynamics.pdf')


def exact_plot(data1, marker="-", label1=""):
    """
    This function plots the matrix elements of the various forms of the dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 

    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0,len(data1[:,:,0]))
            plt.ylim(0,1.1)

            if (i == 2):
                plt.xlabel("t (step)",fontsize=15)
                plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray) / 2)))
                if (j == 0):
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray) / 2)))
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)

            else:
                if (j==0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else: 
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], marker, c='black', linewidth=3.0, markersize=5, markevery=500, label=label1)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)

    plt.savefig('exact.pdf')
    plt.clf()


def plot_exact_NL(data1, data2, color="red", marker="--", marker2="--", label1="", label2=""):
    """
    This function plots the matrix elements of the various forms of the dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 

    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0,len(data1[:,:,0]))
            plt.ylim(0,1.1)

            if (i == 2):
                plt.xlabel("t (step)",fontsize=15)
                plt.xticks(np.arange(min(timearray), max(timearray)+1, 5000))
                if (j == 0):
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, 5000))
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)

            else:
                if (j==0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else: 
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=3.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)

    plt.savefig('exact_and_NL.pdf')
    plt.show()

def plot_NL_L(data1, data2, data3, color="red", color2="blue", marker="--", marker2="-", label1="", label2="", label3=""):
    """
    This function plots the matrix elements of the various forms of the dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 

    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)

    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0,len(data1[:,:,0]))
            plt.ylim(0,1.1)

            if (i == 2):
                plt.xlabel("t (step)",fontsize=15)
                plt.xticks(np.arange(min(timearray), max(timearray)+1, 5000))
                if (j == 0):
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.plot(timearray, data3[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label3)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, 5000))
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else:
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.plot(timearray, data3[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label3)

            else:
                if (j==0):
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.plot(timearray, data3[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label3)
                    plt.yticks(np.arange(0.0, 1.1, 0.5))

                else: 
                    plt.gca().yaxis.set_major_locator(plt.NullLocator())
                    plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i, j, :], "o", c='black', linewidth=3.0, markersize=5, markevery=500, label=label1)
                    plt.plot(timearray, data2[i, j, :], marker, c=color, linewidth=3.0, markersize=10, label=label2)
                    plt.plot(timearray, data3[i, j, :], marker2, c=color2, linewidth=3.0, markersize=10, label=label3)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)

    plt.savefig('exact_NL_L.pdf')
    plt.show() 


def plot_RMSEs(x_axis1, x_axis2, x_axis3, y_axis1, y_axis2,y_axis3, x_title = "time", marker = "o", ub = 2000.0):
    width, length = 8, 5
    plt.figure(0, figsize=(width, length))
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 
    plt.xticks(np.arange(0.0, int(ub)+1, 500)) 
    plt.plot(x_axis1,y_axis1*0.0, marker,c='black',markersize=4.0,linewidth=3.0)
    plt.plot(x_axis1,y_axis1, marker,c='firebrick',markersize=4.0,linewidth=3.0)
    plt.plot(x_axis2,y_axis2,marker,c='royalblue',markersize=3.0,linewidth=3.0)
    plt.plot(x_axis3,y_axis3,marker,c='forestgreen',markersize=3.0,linewidth=3.0)
    plt.xlim(-100.0,ub)
    plt.xlabel(x_title)
    plt.ylabel("RMSE")
    plt.savefig('rmse_analysis.pdf')
    plt.show()
    

def implied_time_scale(time, data, limit, save = True):
    """
    This function plots the implied time scale as a function of time

    Inputs:
    1. time  - time axis
    2. data  - a multidimensional array keeping track of the implied time-scale data
    3. limit - how far to go in the plot; for the 6-state kinetic model this is about 500 time steps.
    4. save  - an option that lets the user save this as a PDF

    Outputs:
    1. a plot of the implied time scales as a function of time. 
    """
    plt.plot(time,data[0,0,:],c='black',label="ITS2")
    plt.plot(time,data[0,1,:],c='firebrick',label="ITS1")
    plt.xlim(0.0,limit)
    plt.xlabel("$t$ (step)")
    plt.ylabel("Implied time scale")
    plt.legend(loc = "center right")
    plt.savefig("ITS_plot.pdf")
    plt.show()


def plot_mem(data1, color="aqua", marker="--", label1="", num=5.0):
    """
    This function plots the reduced transition probability matrix dynamics as a function of time
    
    Inputs:
    1. data - the reduced space transition probability matrix (dim_red x dim_red x tsteps tensor)
    
    Outputs:
    1. n^2 plots; one for the dynamics of each of the matrix elements of the memory kernel
    """
    plt.rc('xtick', labelsize=15) 
    plt.rc('ytick', labelsize=15) 
    
    # number of rows (columns)
    num_rows = data1.shape[0]
    tsteps = data1.shape[-1]
    timearray = np.arange(0, tsteps, 1)
    y0 = np.zeros_like(timearray)

    # set length & width
    width = 12
    length = 12

    # define figure dimensions
    plt.figure(0, figsize=(width, length))

    # plot the matrix elements of the memory kernel
    for i in range(num_rows):
        for j in range(num_rows):

            # useful variables
            subplot_number = num_rows * i + j + 1
            label_index = str(i + 1) + str(j + 1)
            #label_handle1 = r'$R_{' + label_index + '}$(t)'
            #label_handle2 = r'$D_{' + label_index + '}$(t)'

            # create subplots
            plt.subplot(num_rows,num_rows,subplot_number)
            plt.xlim(0, len(data1[0,0,:]))

            if (i == 2):
                if (j == 0):
                    plt.plot(timearray, data1[i,j,:], marker, c=color, linewidth=2.0, markersize=5, markevery=500, label=label1)
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.xticks(np.arange(min(timearray), max(timearray)+1, int(max(timearray) / 2)))

                else:
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i,j,:], marker, c=color, linewidth=2.0, markersize=5, markevery=500, label=label1)
                
            else:
                if (j==0):
                    #plt.gca().xaxis.set_major_locator(plt.NullLocator())
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i,j,:], marker, c=color, linewidth=2.0, markersize=5, markevery=500, label=label1)
                    
                else:
                    plt.xlim(0, len(data1[0,0,:]))
                    plt.plot(timearray, data1[i,j,:], marker, c=color, linewidth=3.0, markersize=5, markevery=500, label=label1)

            if ((i == 2) and (j == 2)):
                plt.legend(loc='upper right', frameon=False)

    plt.tight_layout(pad=num)
    plt.savefig('mem.pdf')
    plt.show()    
