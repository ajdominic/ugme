import numpy as np


def get_U(data):
    """
    This function gets the U matrices needed to propagate the dynamics

    Inputs:
    1. data - the molecular dynamics data 4x4xNt

    Outputs:
    1. U - an array of U matrices that propagae the dynamics 
    """
    # initialize arrays
    U = np.zeros_like(data, float)
    dynamics = np.zeros_like(data, float)

    # get the U matrix at t=0
    U[:, :, 0] = data[:, :, 0]
    #p_ITS = np.zeros([1, int(rd - 1), 21], float)

    # get the U matrices U(n, n-1) = C(n)[C(n-1)]^{-1}
    for k in range(1, data.shape[-1]):
        U[:, :, k] = np.dot(data[:, :, k], np.linalg.inv(data[:, :, k-1]))
        #p_ITS[:, :, k] = msm.psuedo_ITS(temp[:, :, k], rd)

    # reconstruct dynamics to be sure that we did everything correctly
    dynamics[:, :, 0] = U[:, :, 0]

    for k in range(1, dynamics.shape[-1]):
        dynamics[:, :, k] = np.dot(U[:, :, k], dynamics[:, :, k-1])

    # to double check that we have done things correctly
    assert np.allclose(data, dynamics) is True

    return U


def get_Upropagator(data, dim, lb):
    """
    This function builds the propagator needed to calc. long-time dynamics

    Inputs:
    1. data - the time-dependent U_matrices
    2. dim  - the number of rows (columns) of the matrix
    3. lb   - the lower bound for the number of elements used for averaging

    Outputs:
    1. avg  - the average of the data used to calc. long-time dynamics 
    """
    # intialize array
    avg = np.zeros([dim, dim], float)

    # running average
    for k in range(lb, data.shape[-1]):
        avg += data[:, :, k]

    avg /= (data.shape[-1] - lb)

    return avg


def longtime_dynamics(data, U_matrix, file_dt, dim, max_time):
    """
    This function calculates the long-time dynamics.

    Inputs:
    1. data       - the time-dependent data being propagated (size of MD data)
    2. U_matrices - propagates the dynamics from time t-1 to time t
    3. file_dt    - time spacing between data in the provided file
    4. dim        - the number of rows (columns) of the matrix
    5. lb         - the lower bound for the number of elements used for averaging
    6. max_time   - the time that you wish to propagate to
    7. function   - the function used to create the propagator (get_propagator)

    Outputs:
    1. avg  - the average of the data used to calc. long-time dynamics 
    """
    # get the propagator
    #prop = function(U_matrices, dim, lb)

    # initialize the arrays 
    t_axis = np.arange(0.0, max_time+file_dt, file_dt)
    long_dynamics = np.zeros([dim, dim, t_axis.shape[-1]], float)
    
    #ub = ub / file_dt
    for k in range(t_axis.shape[-1]):
        if (k < data.shape[-1]):
        #if (k < ub):
            long_dynamics[:, :, k] = data[:, :, k]

        # if we wish to go beyond the length of MD data provided, we use the U matrices
        else:
            long_dynamics[:, :, k] = np.dot(U_matrix, long_dynamics[:, :, k-1])

    return long_dynamics


def get_dynamics(data, lb, ub, max_time, file_dt, dim, function):
    """
    This function calculates the long-time dynamics.

    Inputs:
    1. data       - the time-dependent data being propagated (size of MD data)
    2. file_dt    - time spacing between data in the provided file
    3. dim        - the number of rows (columns) of the matrix
    4. lb         - the lower bound for the number of elements used for averaging
    5. max_time   - the time that you wish to propagate to
    6. function   - the function used to create the propagator (get_propagator)

    Outputs:
    1. long_dynamics  - the long-time dynamics using the U matrices
    """
    U_matrices = get_U(data)
    long_dynamics = longtime_dynamics(data, U_matrices, file_dt, dim, lb, ub, max_time, function)
    return long_dynamics


def mfpt_calc(M, tau_time):
    """
    This function computes the MFPT at a particular time.
    Inputs:
    1. M        - the TPM at a particular time
    2. tau_time - the time that the TPM is given at
    """
    dim = M.shape[0]

    mfpt = np.zeros_like(M, float)

    # vector of all ones
    vec = np.ones([dim-1, 1], float)
    I = np.identity(dim-1, float)

    for i in range(dim):
        # deletes a row
        temp1 = np.delete(M, (i), axis=0)

        # deletes a column
        temp2 = np.delete(temp1, (i), axis=1)

        # transpose
        #M_tilde = temp2.T
        M_tilde = temp2

        # calculate the prefactor
        pref = np.linalg.inv(I - M_tilde)

        # calculate the mfpt
        calc = tau_time * np.dot(pref, vec)

        # add these back to the mfpt array
        ctr = 0
        for k in range(dim):
            # avoid diagonals
            if (i != k):
                mfpt[i][k] = calc[ctr]
                ctr += 1

    return mfpt
