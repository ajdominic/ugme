import numpy as np
import scipy as sp


def running_ave(matrix, val):
    """
    This function computes the running average of a matrix.
    Inputs:
    1. matrix - the time-dependent matrix we wish to find
                the running ave. of
    2. val    - the time index we wish to start averaging at.

    Outputs:
    1. ave_matrix - the time-dependent matrix output
    """
    dim = matrix.shape[0]
    ave_matrix = np.zeros_like(matrix)
    total = np.zeros([dim, dim])
    for k in range(0, matrix.shape[-1]-1):
        if k > val:
            total += matrix[:, :, k]
            ave_matrix[:, :, k] = total / (k - val)
        else:
            ave_matrix[:, :, k] = matrix[:, :, k]

    return ave_matrix


def num_tcl_gme(exact, dt, tr, long_tsteps, use_exact=False, ave=False, val=20):
    """
    This function computes the time local GME using numerical
    time derivatives of the "MD" data.
    Inputs:
    1. exact       - dynamics in the reduced space (either DDR or MD)
    2. dt          - the time-step corresponding to the exact dynamics
    3. tr          - the generator cut-off time
    4. long_tsteps - the number of time-steps desired to compute

    Outputs:
    1. dictionary containing the generator and the gme data
    """
    # get the R matrix
    R_matrix = np.zeros_like(exact)
    cdot = np.gradient(exact, dt, axis=-1)
    rd, tsteps = R_matrix[0].shape

    # solve for the R_matrix over the time range given
    for k in range(tsteps):
        c_inv = np.linalg.inv(exact[:, :, k])
        R_matrix[:, :, k] = np.dot(cdot[:, :, k], c_inv)

    if ave is True:
        R_matrix = np.copy(running_ave(R_matrix, val))

    # this should recover the exact dynamics
    num_tcl = TCL_evol(R_matrix, rd, tsteps, dt)

    if use_exact is True:
           # next use the R_matrix to do the work for me
        long_num_tcl = TCL_dynamics_with_lt(exact, R_matrix, tr,
                                            dt, long_tsteps)

    else:
        # next use the R_matrix to do the work for me
        long_num_tcl = TCL_dynamics_with_lt(num_tcl, R_matrix, tr,
                                            dt, long_tsteps)

    # build a corresponding time axis
    time_tcl = np.arange(0.0, long_num_tcl.shape[-1], 1) * dt

    # create dictionary
    tcl_dict = {"R":R_matrix, "tcl_gme":num_tcl,
                "long_time":time_tcl, "extended_tcl":long_num_tcl}

    return tcl_dict


def RHS(tn, C_tn, R):
    y = np.dot(R[:, :, tn], C_tn)
    return y


def TCL_evol(generator, rd, tsteps, dt):
    tcl_dyn = np.zeros_like(generator, float)
    tcl_dyn[:, :, 0] = np.identity(rd, float)

    # get C(1)
    k1 = dt * RHS(0, tcl_dyn[:, :, 0], generator)
    dyn_tilde = tcl_dyn[:, :, 0] + k1
    k2 = dt * RHS(1, dyn_tilde, generator)
    tcl_dyn[:, :, 1] = tcl_dyn[:, :, 0] + 0.5 * (k1 + k2)

    for k in range(1, tsteps-1):
        k1 = dt * RHS(k, tcl_dyn[:, :, k], generator)
        dyn_tilde = tcl_dyn[:, :, k] + k1
        k2 = dt * RHS(k+1, dyn_tilde, generator)
        tcl_dyn[:, :, k+1] = tcl_dyn[:, :, k] + 0.5 * (k1 + k2)

    return tcl_dyn


def TCL_dynamics_with_lt(TCL_dynamics, R_matrix, tr_index, dt, tsteps):
    """
    This function propagates the time-local GME using the propagator
    after the generator lifetime has been achieved.

    Inputs:
    1. TCL_dynamics  - the previously calculated TCL dynamics
    2. prop          - the time-local propagator
    3. t_R           - the generator lifetime
    4. tsteps        - the number of steps

    Outputs:
    1. TCL_gme - the time-local GME using the propagator
    """
    if tr_index >= R_matrix.shape[-1]:
        R_infty = R_matrix[:, :, -1]
        tr_index = R_matrix.shape[-1]

    else:
        R_infty = R_matrix[:, :, tr_index]

    prop = get_prop(R_infty, dt)
    rd = R_infty.shape[0]
    #TCL_gme = np.zeros_like(TCL_dynamics, float)
    TCL_gme = np.zeros([rd, rd, tsteps], float)
    for k in range(tsteps):
        if (k < tr_index):
            TCL_gme[:, :, k] = TCL_dynamics[:, :, k]
        else:
            TCL_gme[:, :, k] = np.dot(prop, TCL_gme[:, :, k - 1])

    return TCL_gme


def get_prop(R_infty, dt):
    eig, G = np.linalg.eig(R_infty)
    inv_G = np.linalg.inv(G)
    # abs_eig = np.abs(eig)
    exp_eig = np.exp(eig * dt)
    diag = np.diag(exp_eig)
    # diag = np.diag(abs_eig) * dt
    prop = np.dot(G, np.dot(diag, inv_G))
    return prop


def num_test(sr, timearray, R_matrix):
    """
    This function demonstrates the numerical test
    R(sr) + M / t should plateau like the lagtime plot.

    Inputs:
    1. sr        - the generator cut-off time
    2. timearray - an array keeping track of time
    3. R_matrix  - the generator

    Outputs:
    1. test - an array containing the results that will be plotted.
    """
    dt = 1.0
    R_inf = R_matrix[:, :, sr]
    rd = R_matrix.shape[0]

    # compute the integral
    integral = np.empty([rd, rd], float)
    integral += (R_matrix[:, :, 0] + R_matrix[:, :, sr]) / 2
    for k in range(1, sr):
        integral += R_matrix[:, :, k]
    integral *= dt

    tsteps = len(timearray)
    test = np.empty_like(R_matrix, float)

    # the numerical test for R(sr)
    RITS = np.empty([1, int(rd - 1), tsteps], dtype='float')
    for k in range(1, tsteps):
        term1 = R_inf * (1 - sr / k)
        term2 = integral / k
        test[:, :, k] = 1 / (term1 + term2)
        RITS[:, :, k] = RITS_test(test[:, :, k], rd)

    return test, RITS


def RITS_test(tau, dim_red):
    """
    This function returns ITS(t) for each of N = (dimension - 1) eigenvalues

    Inputs:
    1. tau     - the TPM in the reduced space (matrix)
    2. time    - the time-step of this particular propagation (scalar)
    3. dim_red - the dimension of the reduced space (scalar)

    Outputs:
    1. ITS - 1 x num array containing the implied time scale value at this time
             step for the N-1 eigenvalues (1 x num dimensional array)
    """
    # the number of eigenvalues not equal to 1
    num = dim_red - 1

    # get the eigenvalues of tau and sort them
    eig, unitary = np.linalg.eig(tau)
    eig = np.sort(eig, axis=None)

    # loop to compute and store the ITS values
    ITS = np.zeros(num)
    for j in range(0, num):
        ITS[j] = 1 / eig[j]

    return ITS
        