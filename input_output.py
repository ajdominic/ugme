import numpy as np
import pickle


def load_data(filename):
    """Loads data from a pickle file"""

    with open(filename, 'rb') as f:
        data = pickle.load(f)

    return data


def save_data(data, filename):
    """Saves data to a picklefile."""

    with open(filename, 'wb') as f:
        pickle.dump(data, f)


def get_count_matrix(filename, n=0):
    """
    This function produces the count matrix based on data given.
    Inputs:
    1. filename - full path to the transition count matrix data

    Outputs:
    1. TCM - the transition count matrix
    """
    # open the file and store it
    file = open(filename, "r")
    line = file.readlines()[n]
    file.close()

    # store the line as a list and obtain dimensions
    temp_list = line.split()
    size = int(len(temp_list))
    dim = int(size**0.5)

    # create the transition count matrix
    TCM = np.zeros(size, float)
    for k in range(size):
        TCM[k] = float(temp_list[k])

    # reshape the matrix from (1, size) to (dim, dim)
    #TCM = TCM[1:].reshape((dim, dim))
    TCM = TCM.reshape((dim, dim))

    return TCM


def new_assignment_matrix(msm_dim):
    """Code to randomly construct the assignment matrix."""
    pre_assnt = np.zeros([msm_dim, 2], int)
    for i in range(msm_dim):
        pre_assnt[i, 1] = int(np.floor(np.random.rand() * 4.0))
        pre_assnt[i, 0] = i
    return pre_assnt


def get_assnt_matrix(filename):
    """
    This function produces the assignment matrix based on data given.
    Inputs:
    1. filename - full path to the transition count matrix data

    Outputs:
    1. assnt - the assignment matrix
    """
    # open the file and store it
    file = open(filename, "r")
    temp = file.read()
    file.close()

    # store as a list and get dimensions
    raw_text = temp.split()
    size = int(len(raw_text))
    dim = int(size / 2)

    # copy the list to a numpy array
    pre_assnt = np.zeros(size, int)
    for k in range(size):
        pre_assnt[k] = int(raw_text[k])

    # reshape for aesthetic purposes
    pre_assnt = pre_assnt.reshape((dim, 2))
    #pre_assnt = new_assignment_matrix(dim)

    # only keeping column 2; this has all the necessary info
    assnt_values = np.zeros(dim, int)
    for k in range(dim):
        assnt_values[k] = pre_assnt[k][1]

    # the dimension of the reduced space is the max value in col 2
    rd = max(assnt_values) + 1

    # create the skeleton for the assignment matrix
    assnt = np.zeros([dim, rd], float)
    for i in range(dim):
        # A_iJ = 1 if i \in J and A_iJ = 0 if i \notin J
        J = assnt_values[i]
        assnt[i][J] = 1.0

    return assnt


def load_md(file_name, md_fn, er_fn, rd=4, init=True):
    """
    This function loads md data with error bars and saves to a pickle file.

    Inputs:
    1. file_name - the name of the file you wish to read in data from
    2. md_fn     - the name of the file you wish to save the md data to 
    3. er_fn     - the name of the file you wish to save the errors to
    4. rd        - the dimension of the space that you are projecting onto

    Outputs:
    1. the md data is saved to a pickle file
    2. the error data is saved to a pickle file
    """
    # get number of lines and initialize arrays
    num_lines = sum(1 for line in open(file_name))
    #md_data = np.zeros([rd, rd, num_lines+1], float)
    #er_bars = np.zeros([rd, rd, num_lines+1], float)

    if init is True:
        p = 1
        md_data = np.zeros([rd, rd, num_lines+p], float)
        er_bars = np.zeros([rd, rd, num_lines+p], float)
        md_data[:, :, 0] = np.identity(rd)

    else:
        p=0
        md_data = np.zeros([rd, rd, num_lines+p], float)
        er_bars = np.zeros([rd, rd, num_lines+p], float)

    # useful indexing parameters
    rd2 = rd**2 + 1
    rd3 = 2 * rd2 - 1

    # reading in file
    file = open(file_name, "r")
    lines = file.readlines()
    time = []
    for k in range(num_lines):
        # read in line
        line = lines[k]
        temp_list = line.split()
        temp_array = np.array(temp_list, float)
        time.append(temp_array[0])

        temp2 = temp_array[1:rd2]
        temp3 = temp_array[rd2:rd3]
        #temp3 = temp_array[1+rd2:]
        #temp3 = temp_array[1+rd2:rd3]
        #temp2 = temp_array[1:17]
        #temp3 = temp_array[18:34]

        # store the md data
        md = temp2.reshape((rd, rd))
        md_data[:, :, k+p] = md

        # store the md error bars
        er_temp = temp3.reshape((rd, rd))
        er_bars[:, :, k+p] = er_temp

    # close the file
    file.close()

    # save the data
    save_data(md_data, md_fn)
    save_data(er_bars, er_fn)

    # return the data
    return md_data, er_bars

def load_md_no_er(file_name, md_fn, rd=4, init=True):
    """
    This function loads md data with error bars and saves to a pickle file.

    Inputs:
    1. file_name - the name of the file you wish to read in data from
    2. md_fn     - the name of the file you wish to save the md data to 
    3. er_fn     - the name of the file you wish to save the errors to
    4. rd        - the dimension of the space that you are projecting onto

    Outputs:
    1. the md data is saved to a pickle file
    2. the error data is saved to a pickle file
    """
    # get number of lines and initialize arrays
    num_lines = sum(1 for line in open(file_name))
    #md_data = np.zeros([rd, rd, num_lines+1], float)
    #er_bars = np.zeros([rd, rd, num_lines+1], float)

    if init is True:
        p = 1
        md_data = np.zeros([rd, rd, num_lines+p], float)
        md_data[:, :, 0] = np.identity(rd)

    else:
        p=0
        md_data = np.zeros([rd, rd, num_lines+p], float)

    # useful indexing parameters
    rd2 = rd**2 + 1
    rd3 = 2 * rd2 - 1

    # reading in file
    file = open(file_name, "r")
    lines = file.readlines()
    time = []
    for k in range(num_lines):
        # read in line
        line = lines[k]
        temp_list = line.split()
        temp_array = np.array(temp_list, float)
        time.append(temp_array[0])

        if rd==4:
            temp2 = temp_array[0:rd2]
        else:
            temp2 = temp_array[1:rd2]
        #temp3 = temp_array[rd2:rd3]
        #temp3 = temp_array[1+rd2:]
        #temp3 = temp_array[1+rd2:rd3]
        #temp2 = temp_array[1:17]

        # store the md data
        md = temp2.reshape((rd, rd))
        md_data[:, :, k+p] = md


    # close the file
    file.close()

    # save the data
    save_data(md_data, md_fn)

    # return the data
    return md_data