import numpy as np


def rmse(array1, array2):
    """
    Inputs:
    1. array1  - true GME dynamics (from direct dim. red.)
    2. array2  - GME dynamics at given kernel (or generator) lifetime

    Outputs:
    1. sqrt(MSE) - the root mean square error at a particular lifetime
    """
    N = array1.shape[-1]
    dim = array1.shape[0]
    total = 0.0

    for k in range(N):
        difference = array1[:, :, k] - array2[:, :, k]
        difference_squared = difference**2
        total += np.sum(difference_squared)

    total = (total / (N * dim**2))**0.5

    return total


def TC_rmse(integrator, timearray, mem, exact, omega, tsteps, dt, alpha=1.0):
    TC_rmse_array = np.zeros_like(timearray, float)

    for k in range(len(timearray)):
        conv_trial = integrator(mem, exact[:, :, 0], omega,
                                tsteps, k, alpha, dt)
        rmse_val = rmse(exact, conv_trial)
        TC_rmse_array[k] = rmse_val

    return TC_rmse_array


def TCL_rmse(integrator, timearray, exact, local_gme, R_matrix, dt, tsteps):
    TCL_rmse_array = np.zeros_like(timearray, float)

    for k in range(len(timearray)):  
        TCL_temp = integrator(local_gme, R_matrix, k, dt, tsteps)
        rmse_val_R = rmse(exact, TCL_temp)
        TCL_rmse_array[k] = rmse_val_R

    return TCL_rmse_array


def RDMSM_rmse(integrator, timearray, exact, t_max, dt):
    rdmsm_rmse_array = np.zeros_like(timearray, float)
    
    #time = np.arange(0, t_max+1, 5)
    for k in range(1, int(len(timearray))-1):
    #for k in range(1, int(len(time))):
        # in this loop, each k is a reduced space lagtime
        #lt_index = int(k/dt) + 1
        trial_lt = k
        t_temp, rdmsm_temp = integrator(trial_lt, exact, t_max, rd, dt)
        rmse_val = msm_rmse(exact, rdmsm_temp, trial_lt)
        rdmsm_rmse_array[k] = rmse_val

    return rdmsm_rmse_array


def msm_rmse(exact_data, msm_data, lt):
    """
    This function computes the RMSE for a MSM with respect
    to the exact dynamics. This is a separate function because
    the two arrays differ in size.
    Inputs:
    1. array1 - Exact, dimensionality reduced dynamics
    2. array2 - The reduced MSM dynamics
    3. lt     - The lag time that the RDMSM was computed at
    4. tsteps - The number of steps
    """
    N = msm_data.shape[-1]
    total = 0.0

    for k in range(2, N):
        j = int(lt * k)
        if (j < exact_data.shape[-1] and j < N):
            difference = exact_data[:, :, j] - msm_data[:, :, k]
            difference_squared = difference**2
            total += np.sum(difference_squared)

    total = (total / N)**0.5

    return total

def threshold(error, RMSE):
    """
    Returns the index of the first instance of the RMSE being below a given
    error threshold.

    Inputs:
    1. error - a scalar value indicating the threshold of error
    2. RMSE  - RMSE is a 1xtsteps array that contains the error of GME
               dynamics compared to the exact dynamics as a function of time
    """
    ctr = 0
    value = False
    if RMSE[0] < error:
        value = True

    while ((value is False) and (ctr < RMSE.shape[0])):
        ctr += 1
        if RMSE[ctr] < error:
            value = True

    return ctr, RMSE[ctr]
    

#TCL_rmse_array = np.zeros_like(timearray, float)
#for k in range(len(timearray)):
    #R_infty = R_matrix[:, :, tr_index]
    #prop_temp = get_prop(R_infty, dt)
    #TCL_temp = tcl.TCL_dynamics_with_lt(exact, prop_temp, tr_index, tsteps)
    #rmse_val_R = msm.rmse(exact, TCL_temp, dt)
    #TCL_temp = tcl.TCL_dynamics_with_lt(local_gme, R_matrix, k, dt, tsteps)
    #rmse_val_R = er.rmse(exact, TCL_temp)
    #TCL_rmse_array[k] = rmse_val_R
