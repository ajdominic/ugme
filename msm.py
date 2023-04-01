import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
import pickle


def markovian_prop(prop, tsteps, rd):
    """
    This function returns the dynamics corresponding to a given 'Markovian'
    propagator. I generally do not use this to calculate MSMs
    Input:
    1. prop - a [rd, rd] matrix that is of the form e^(Ω+K(0)) or e^(Ω+intK)}
    2. tsteps - the number of timesteps
    3. rd - the reduced space dimension (usually 3 or 4)

    Output:
    1. dynamics - Markovian dynamics corresponding to the prop [rd, rd, tsteps]
    """
    # initialize shiz
    dynamics = np.zeros([rd, rd, tsteps], float)
    dynamics[:, :, 0] = np.eye(rd)

    # do the evolution
    for k in range(1, tsteps):
        dynamics[:, :, k] = np.dot(prop, dynamics[:, :, k-1])

    return dynamics


def msm_folding_times(threshold, exact, tmax, rd=3, dt=1, exiting=1, to=2):

    folding_time_cut_offs = []
    folding_time_list = []
    ctrs = []

    sample_lts = np.arange(2, exact.shape[-1], 1)

    for lt in sample_lts:
        # compute dynamics with the U matrix; EX: (exact, ave2(U_matrix, 20, 5), file_dt, rd, 2000, 25)
        steps = int(tmax / lt)
        msm_temp = msm_calc1(exact, lt, dt, steps)

        # compute corresponding MFPT
        Umfpt = np.zeros_like(msm_temp)
        for k in range(1, Umfpt.shape[-1]):
            Umfpt[:, :, k] = m.mfpt_calc(msm_temp[:, :, k], 1)

        M = np.copy(Umfpt[to, exiting, :])
        last = M[-1]
        diff = np.abs(M - last)

        ctr = 0
        test = False
        while (test is False) and (ctr * lt < exact.shape[-1]):
        #while ctr < tmax:
        #for ctr in range(tmax):
            if diff[ctr] >= threshold:
                ctr += 1
            elif diff[ctr] < threshold:
                folding_time_cut_offs.append(lt)
                folding_time_list.append(M[ctr]*ctr*lt/1000)
                ctrs.append(ctr)
                test = True

        if ctr * lt >= exact.shape[-1]:
            folding_time_cut_offs.append(lt)
            folding_time_list.append(M[ctr]*ctr*lt/1000)
            ctrs.append(ctr)

    return folding_time_cut_offs, folding_time_list, ctrs


def folding_times(U_matrix, threshold, exact, tmax, rd=3, dt=1, tr=15, exiting=1, to=2):

    folding_time_cut_offs = []
    folding_time_list = []
    ctrs = []

    sample_tlrs = np.arange(tr+1, exact.shape[-1], 1)
    for tlr in sample_tlrs:
        # compute dynamics with the U matrix; EX: (exact, ave2(U_matrix, 20, 5), file_dt, rd, 2000, 25)
        U_temp = longtime_dynamics(exact, ave2(U_matrix, tr, tlr-tr), dt, rd, tmax, tlr)

        # compute corresponding MFPT
        Umfpt = np.zeros_like(U_temp)
        for k in range(1, Umfpt.shape[-1]):
            Umfpt[:, :, k] = m.mfpt_calc(U_temp[:, :, k], 1)

        M = np.copy(Umfpt[to, exiting, :])
        last = np.copy(M[-1])
        diff = np.abs(M - last)

        ctr = 0
        test = False
        while (test is False) and (ctr < tmax):
        #while ctr < tmax:
        #for ctr in range(tmax):
            if diff[ctr] >= threshold:
                ctr += 1
            elif diff[ctr] < threshold:
                folding_time_cut_offs.append(tlr)
                folding_time_list.append(M[ctr]*ctr/1000)
                ctrs.append(ctr)
                test = True

    return folding_time_cut_offs, folding_time_list, ctrs


def ave(U_matrix, val, num):
    """
    returns the average U matrix
    Inputs:
    1. U_matrix - the full time-dependent U matrix
    2. val - the index to begin averaging at
    3. num - end of averaging window
    """
    return np.average(U_matrix[:, :, val:val+num], axis=-1)


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
    ave_matrix = np.zeros_like(matrix)
    average = ave2_v2(matrix, val)
    ave_matrix[:, :, 0] = matrix[:, :, 0]
    for k in range(1, matrix.shape[-1]-1):
        if k > val:
            ave_matrix[:, :, k] = np.copy(average)
        else:
            ave_matrix[:, :, k] = matrix[:, :, k]

    return ave_matrix


def longtime_dynamics(data, U_matrix, file_dt, dim, tsteps, tr):
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
    t_axis = np.arange(0.0, tsteps, 1.0) * file_dt
    long_dynamics = np.zeros([dim, dim, tsteps], float)
    
    #ub = ub / file_dt
    for k in range(tsteps):
        if (k < tr):
        #if (k < ub):
            long_dynamics[:, :, k] = data[:, :, k]

        # if we wish to go beyond the length of MD data provided, we use the U matrices
        else:
            long_dynamics[:, :, k] = np.dot(U_matrix, long_dynamics[:, :, k-1])

    return long_dynamics


def trap_int(matrix, ub_index, dt):
    """
    Integrates dynamics according to a trapzoidal rule
    Inputs:
    1. matrix - dynamical object [rd, rd, tsteps]
    2. ub_index - upper limit of the integral
    3. dt - timestep

    Outputs:
    1. dt * integral - the integrated form of original matrix
    """
    dim = matrix.shape[0]
    integral = np.zeros([dim, dim])

    # handle the end points
    integral += (matrix[:, :, 0] + matrix[:, :, -1]) / 2.0

    # handle the others
    for k in range(1, ub_index):
        integral += matrix[:, :, k]

    # make sure everyone gets a dt
    return dt * integral


def mik_calc(matrix): # put this into the msm code
    """
    mik := mean integral of the memory kernel

    Inputs:
    1. matrix - memory kernel or time-local generator
                with shape [d, d, tsteps]

    Outputs:
    1. mik - the mik calculation
    """
    rd = matrix.shape[0]
    temp = np.zeros_like(matrix)
    for k in range(1, temp.shape[-1]):
        val = trap_int(matrix, k, 1)
        temp[:, :, k] = val**2 # squares all matrix elems

    mik = np.zeros(temp.shape[-1])
    for k in range(1, temp.shape[-1]):
        mik[k] = np.sqrt(np.sum(temp[:, :, k])) / rd

    return mik


def num_tc_gme(exact, dt, tk, t_max, alpha=1.0):
    """
    This function computes the time convolution GME using numerical
    time derivatives of the "MD" data.
    Inputs:
    1. exact - dynamics in the reduced space (either DDR or MD)
    2. dt    - the time-step corresponding to the exact dynamics
    3. tk    - the kernel cut-off time

    Outputs:
    1. dictionary containing the kernels and the gme data
    """
    # numerical stuff
    cdot = np.gradient(exact, dt, axis=-1)
    cddot = np.gradient(cdot, dt, axis=-1)

    # some parameters
    rd, tsteps = cdot[0].shape

    # get the time derivative of C at zero (omega)
    c0 = exact[:, :, 0]
    cd0 = cdot[:, :, 0]

    # initialize the aux
    auxk1 = np.zeros_like(cdot)
    auxk3 = np.zeros_like(cdot)

    for k in range(cdot.shape[-1]):
        # first build k1
        term1 = cddot[:, :, k]
        term2 = - np.dot(cdot[:, :, k], cd0) - np.dot(cd0, cdot[:, :, k])
        term3 = np.dot(np.dot(cd0, exact[:, :, k]), cd0)

        # add to k1 and k3
        auxk1[:, :, k] = term1 + term2 + term3
        auxk3[:, :, k] = np.dot(cd0, exact[:, :, k]) - cdot[:, :, k]

    # get the memory kernel
    kernel = gen_mem(auxk1.shape[-1], dt, rd, auxk1, auxk3)
    num_tc = gme_evolution(kernel, c0, cd0, t_max, tk, alpha, dt)

    tc_dict = {"k1":auxk1, "k3":auxk3, "kernel":kernel, "tc_gme":num_tc}

    return tc_dict


def get_EXACT(TCM_file, assnt, num_lines, file_dt, sym=True):
    """
    This function gets the exact dynamics from direct dimensionality reduction.

    Inputs:
    1. TCM_file  - the name of the file containing the TCM data
    2. assnt     - the assignment matrix used to create the projector
    3. num_lines - number of lines in the TCM file
    4. rd        - the dimension of the space being projected onto
    5. file_dt   - dt used in the file (resolution of the data provided) 
    6. sym       - used if the TCM data has not been symmetrized

    Outputs:
    1. t_EXACT - time array corresponding to the dimensionality reduced data
    2. EXACT   - the exact dynamics from direct dimensionaliry reduction
    """
    # get dimensions of the assignment matrix
    dim, rd = assnt.shape

    # initialize arrays
    t_EXACT = [0.0]
    EXACT = np.zeros([rd, rd, num_lines], float)
    EXACT[:, :, 0] = np.identity(rd, float)

    # get the PAD component of the projector
    C_temp = get_count_matrix(TCM_file, num_lines-1, True)
    TPM_temp = normalize_columns(C_temp)
    pop_eq = propagate(TPM_temp, num_lines, dim)
    PAD = get_projector(assnt, pop_eq)

    for k in range(1, num_lines):
        t_EXACT.append(float(k))

        # create TPM
        C = get_count_matrix(TCM_file, k-1, True)

        if sym is True:
            C = (C + C.T) / 2.0

        # construct TPM from TCM
        TPM = normalize_columns(C)

        # get data
        EXACT[:, :, k] = np.dot(assnt.T, np.dot(TPM, PAD))

    t_EXACT = np.array(t_EXACT, float) * file_dt

    return [t_EXACT, EXACT]


def get_aux_kernels(ders, assnt, TPAD):
    """
    This function returns the auxiliary kernels K1(t) and K3(t)

    Inputs:
    1. td0  - the first derivative of tau at t=0 (matrix)
    2. A    - the aggregation matrix (matrix)
    3. M    - the first derivative of T at t=0 (matrix)
    4. M2   - the first derivative of T at t=0 squared (matrix)
    5. TPAD - the product of T(t) * PI * A * D^{-1} (matrix)

    Outputs:
    1. K1 - auxiliary kernel K^(1)(t) (matrix)
    2. K3 - auxiliary kernel K^(3)(t) (matrix)
    """
    M = ders[0]
    M2 = np.dot(M, M)
    omega = ders[1]

    # compute the first and second derivatives of tau(t)
    c_dot = np.dot(assnt.T, np.dot(M, TPAD))
    c_ddot = np.dot(assnt.T, np.dot(M2, TPAD))
    c_t = np.dot(assnt.T, TPAD)

    # collect terms needed to compute the auxiliary kernels
    term_2 = np.dot(c_dot, omega)
    term_3 = np.dot(omega, c_dot)
    temp = np.dot(omega, c_t)
    term_4 = np.dot(temp, omega)

    # term_1 is tau_2nd_deriv
    K1 = c_ddot - (term_2 + term_3) + term_4
    K3 = - c_dot + temp

    return K1, K3


def get_count_matrix(filename, n=0, shift=False):
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
    if shift is True:
        TCM = TCM.reshape((dim, dim))
    else:
        TCM = TCM[1:].reshape((dim, dim))
    
    return TCM


def msm_dyn(lt, EXACT, t_max, rd, dt):
    """
    Computes the msm dynamics in the reduced space

    Inputs:
    1. lt    - the lagtime you wish to use to construct msm
    2. EXACT - the exact dynamics obtained from ddr
    3. t_max - length of time for the calc
    4. rd    - dimension of the space being projected onto
    5. dt    - dt of the EXACT dynamics
    """
    # get time info
    t_axis = np.arange(0.0, t_max+1, lt)
    tsteps = t_axis.shape[-1]

    # propagator based on a given lt
    #lt_index = int(lt / dt)
    lt_index = lt
    prop = EXACT[:, :, lt_index]

    # generate the msm dynamics
    msm_dyn = np.zeros([rd, rd, tsteps], float)
    msm_dyn[:, :, 0] = np.identity(rd, float)
    for k in range(1, tsteps):
        msm_dyn[:, :, k] = np.dot(prop, msm_dyn[:, :, k-1])

    # these need to be equal by construction
    assert np.allclose(msm_dyn[:, :, 1], prop)

    # create an exact array that is parallel to the msm array
    exact_par = np.zeros_like(msm_dyn, float)

    ctr = 0
    for k in range(exact_par.shape[-1]):
        #new_index = k * lt_index
        if (ctr < EXACT.shape[-1]):
            exact_par[:, :, k] = EXACT[:, :, ctr]
            ctr += lt_index-1

    return t_axis, msm_dyn, exact_par


def symmetrize(matrix):
    return (matrix + matrix.T) / 2.0


def normalize_columns(matrix):
    """
    This function normalizes the columns of the given matrix
    by dividing each element of a given column by the sum of
    the elements of that column. Specifying an order of 1 in
    the linalg.norm indicates that we are dividing by the sum
    of the column. Summing down a column is specified by
    setting axis = 0.
    """
    return matrix / np.linalg.norm(matrix, ord=1, axis=0, keepdims=True)


def propagate(T, n, dim):
    """
    Propagates the transition probability matrix (TPM) n times

    Inputs:
    1. T   - the transition probability matrix (TPM)
    2. n   - the number of times that the TPM should be propagated
    3. dim - the number of columns (rows) of the T matrix

    Output:
    1. The propagated TPM
    2. equilibrium population, pi vector
    """
    Tn = np.linalg.matrix_power(T, n)
    Tn1 = np.linalg.matrix_power(T, n + 1)
    test = np.allclose(Tn, Tn1)
    while(test is False):
        n += 1
        Tn = np.linalg.matrix_power(T, n)
        Tn1 = np.linalg.matrix_power(T, n + 10)
        test = np.allclose(Tn, Tn1)
    v = np.zeros(dim)
    v[0] = 1.0

    pop_eq = np.dot(Tn1, v)

    return pop_eq


def discretize(T, dt, ms_lt):
    """
    This function correctly discretizes the TPM.
    For some reason scipy's exp is more stable with smaller timesteps
    Inputs:
    1. T     - The many-state TPM at a given lagtime, ms_lt
    2. dt    - The time-step we wish to discretize to.
    3. ms_lt - The many-state lagtime that parametrizes the current TPM
    """
    eig, P = np.linalg.eig(T)
    P_inv = np.linalg.inv(P)
    diag = np.diag(np.abs(eig))
    exponent = dt / float(ms_lt)
    temp = diag**exponent
    disc_TPM = np.dot(P, np.dot(temp, P_inv))
    return disc_TPM


def get_M(T, dt, ms_lt):
    """
    This function takes in the trasition probability (TPM) matrix at T(1)
    and returns the [dot{T}_{j,k}](0), [<j|-M|k rho>]

    Inputs:
    1. T is the TPM at time t = 1, T_{j,k}(1) = <j|e^{-L}|k rho>pi_k^{-1};
        T(1) = e^{-M} (matrix)

    Outputs:
    1. TPM_d1 - the first derivative of the many-state TPM evaluated at t=0
                normalized by the many-state lagtime chosen.
    2. TPM    - the TPM at the first time step (0 + dt)
    """

    # find the invertible matrix P that diagonalizes T(1), P^{-1}e^{-M}P = D
    eig, P = np.linalg.eig(T)
    P_inv = np.linalg.inv(P)

    # get log of eigenvalues and put on the diagonal of a matrix
    log_evs = np.diag(np.log(np.abs(eig)))
    #log_evs = np.diag(np.log(eig))
    M = - np.dot(P, np.dot(log_evs, P_inv))
    TPM_d1 = - M / ms_lt
    TPM = expm(TPM_d1 * dt)
    #TPM = discretize(T, dt, ms_lt)

    return TPM_d1, TPM


def get_projector(assnt, pop_eq):
    """
    This function computes the PAD part of the projector
    """
    PI = np.diag(pop_eq)
    PI_inv = np.linalg.inv(PI)
    DN = np.dot(assnt.T, np.dot(PI, assnt))
    PAD = np.dot(np.dot(PI, assnt), np.linalg.inv(DN))
    return PAD


def DR_exact_dynamics(TPM, projector, ders, tsteps, dt):
    """
    This function only computes the dimensionality reduced dynamics
    obtained via direct DR.
    """
    # get info for the projector
    assnt = projector[0]
    PAD = projector[1]
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    M = ders[0]
    omega = ders[1]

    dynamics = np.zeros([rd, rd, tsteps], float)
    generator = np.zeros_like(dynamics, float)
    # get exact dynamics, memory kernel, and generator
    Tn_loop = np.identity(dim, float)
    dynamics[:, :, 0] = np.dot(assnt.T, PAD)
    TPAD = np.dot(Tn_loop, PAD)

    # get the generator at t = 0
    generator[:, :, 0] = omega

    for k in range(1, tsteps):
        # propagating the many-state model
        Tn_loop = np.dot(Tn_loop, TPM)

        # dimensionality reduction
        TPAD = np.dot(Tn_loop, PAD)
        dynamics[:, :, k] = np.dot(assnt.T, TPAD)
        inv = np.linalg.inv(dynamics[:, :, k])

        temp = np.dot(M, Tn_loop)
        red_prime = np.dot(assnt.T, np.dot(temp, PAD))
        generator[:, :, k] = np.dot(red_prime, inv)

    return dynamics, generator


def mem_k_dynamics(TPM, projector, ders, tsteps, dt):
    """
    This function only computes the dimensionality reduced dynamics
    obtained via direct DR.
    """
    assnt = projector[0]
    PAD = projector[1]
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    M = ders[0]
    omega = ders[1]

    dynamics = np.zeros([rd, rd, tsteps], float)
    auxk1 = np.zeros_like(dynamics, float)
    auxk3 = np.zeros_like(dynamics, float)

    # get exact dynamics, memory kernel, and generator
    Tn_loop = np.identity(dim, float)
    dynamics[:, :, 0] = np.dot(assnt.T, PAD)
    TPAD = np.dot(Tn_loop, PAD)
    auxk1[:, :, 0], auxk3[:, :, 0] = get_aux_kernels(ders, assnt, TPAD)

    for k in range(1, tsteps):
        # propagating the many-state model
        Tn_loop = np.dot(Tn_loop, TPM)

        # dimensionality reduction
        TPAD = np.dot(Tn_loop, PAD)
        auxk1[:, :, k], auxk3[:, :, k] = get_aux_kernels(ders, assnt, TPAD)
        dynamics[:, :, k] = np.dot(assnt.T, TPAD)

    return dynamics, auxk1, auxk3


def get_long_aux(tcm_file, assnt, time, dt, ms_lt, sym=True):
    timearray = np.arange(0.0, time, dt)
    tsteps = timearray.shape[-1]
    alpha = 1.0
    tk = tsteps
    tk_index = int(tk / dt)
    
    # obtain count matrix and construct projectors and derivates
    C = get_count_matrix(tcm_file, int(ms_lt)-1)

    # symmetrize the count matrix
    if sym is True:
        C = (C + C.T) / 2.0

    TPM = normalize_columns(C)
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    pop_eq = propagate(TPM, int(t_max), dim)
    TPM_d1, TPM = get_M(TPM, dt, ms_lt)
    PAD = get_projector(assnt, pop_eq)
    projector = [assnt, PAD]
    omega = np.dot(np.dot(assnt.T, TPM_d1), PAD)
    ders = [TPM_d1, omega]

    # calculate K1(t), K3(t), and C^EXACT(t)
    exact, k1, k3 = mem_k_dynamics(TPM, projector, ders, tsteps, dt)
    # R_matrix = msm.DR_exact_dynamics(TPM, projector, ders, tsteps, dt)[1]

    """
    plt.figure(0, figsize=(12, 12))
    for k in range(rd):
        for j in range(rd):
            subplot_number = rd * k + j + 1
            plt.subplot(rd, rd, subplot_number)
            plt.plot(timearray, k1[k, j, :], label="k1")
            plt.plot(timearray, k3[k, j, :], label="k3")
            plt.xlim(0.0, time)
        plt.tight_layout(pad=0)
    plt.legend(loc="lower right")
    plt.savefig("./output/mem_kernels.pdf")
    plt.show()"""

    return exact, k1, k3


def get_exact(tcm_file, assnt, time, dt, file_dt, ms_lt, sym=True):
    timearray = np.arange(0.0, time, dt)
    tsteps = timearray.shape[-1]
    t_max = time
    alpha = 1.0
    lt_index = int(ms_lt / file_dt)

    # obtain count matrix and construct projectors and derivates
    C = get_count_matrix(tcm_file, lt_index-1, True) # changed last argument to true 11/16/22

    # symmetrize the count matrix
    if sym is True:
        C = (C + C.T) / 2.0

    TPM = normalize_columns(C)
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    pop_eq = propagate(TPM, int(t_max), dim)
    TPM_d1, TPM = get_M(TPM, dt, ms_lt)
    PAD = get_projector(assnt, pop_eq)
    projector = [assnt, PAD]
    omega = np.dot(np.dot(assnt.T, TPM_d1), PAD)
    ders = [TPM_d1, omega]

    # calculate K1(t), K3(t), and C^EXACT(t)
    exact, k1, k3 = mem_k_dynamics(TPM, projector, ders, tsteps, dt)
    R = DR_exact_dynamics(TPM, projector, ders, tsteps, dt)[1]
    mem = gen_mem(tsteps, dt, rd, k1, k3)

    return exact, k1, k3, mem, R, omega


def get_only_exact(tcm_file, assnt, time, dt, file_dt, ms_lt, sym=True):
    """this function returns the continuous exact dynamics"""
    alpha = 1.0

    # time stuff
    timearray = np.arange(0.0, time, dt)
    tsteps = timearray.shape[-1]

    # kernel cut-off time and MSM lagtime
    tk = tsteps
    tk_index = int(tk / file_dt)
    lt_index = int(ms_lt / file_dt)
    
    # obtain count matrix and construct projectors and derivates
    C = get_count_matrix(tcm_file, lt_index-1)

    # symmetrize the count matrix
    if sym is True:
        C = (C + C.T) / 2.0

    # get TPM and equilibrium populations
    TPM = normalize_columns(C)
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    pop_eq = propagate(TPM, int(time), dim)

    # the 
    TPM_d1, TPM = get_M(TPM, dt, ms_lt)
    PAD = get_projector(assnt, pop_eq)
    projector = [assnt, PAD]
    omega = np.dot(np.dot(assnt.T, TPM_d1), PAD)
    ders = [TPM_d1, omega]

    # calculate K1(t), K3(t), and C^EXACT(t)
    exact = mem_k_dynamics(TPM, projector, ders, tsteps, dt)[0]

    return exact


def RD_MSM(prop, t_max, rdlt):
    """
    This propagates the MSM in the reduced space.

    Inputs:
    1. prop  - the propagator (from the exact dynamics at
              a particular timestep)
    2. t_max - the time for the simulation
    3. rdlt  - the reduced space lagtime (rdlt >= mslt)

    Outputs:
    1. rdmsm - the reduced space msm dynamics
    """
    timearray = np.arange(0.0, t_max+1, rdlt)
    N = len(timearray)
    rd = prop.shape[0]

    # initialize the rdmsm
    rdmsm = np.zeros([rd, rd, N], float)
    rdmsm[:, :, 0] = np.identity(rd, float)

    # simple propagation
    for i in range(1, N):
        rdmsm[:, :, i] = np.dot(prop, rdmsm[:, :, i-1])

    return rdmsm


def lt_plot_analysis(exact, rd, spacing, tsteps):
    """
    old code
    """
    lagtimes = np.arange(spacing, tsteps + 1, spacing)
    rmse_values = np.empty_like(lagtimes, float)

    ctr = 0
    while (ctr < len(lagtimes) + 1):
        lt = lagtimes[ctr]
        prop = exact[:, :, lt]
        temp_msm = dynamics_at_lt(rd, lt, tsteps, prop)
        value = rmse(temp_msm, exact, spacing, tsteps)
        rmse_values[ctr] = value

    # for j, k in zip(range(len(lagtimes)), lagtimes):
        # if ((k != 0) and (k < gme.shape[-1])):
            # temp_msm, temp_plot = dynamics_at_lt(exact, rd, k, tsteps, tau[:, :, k])
            # value = rmse(temp_msm, exact, k)
            # rmse_values[j] = value

    return lagtimes, rmse_values


def exact_and_ITS(tcm_file, assnt, time, dt, ms_lt, timearray, sym=True):
    t_max = time
    tsteps = int((t_max + 1) / dt) + 1
    alpha = 1.0
    tk = tsteps
    tk_index = int(tk / dt)
    
    # obtain count matrix and construct projectors and derivates
    C = get_count_matrix(tcm_file, int(ms_lt)-1)

    # symmetrize the count matrix
    if sym is True:
        C = (C + C.T) / 2.0

    TPM = normalize_columns(C)
    dim = assnt.shape[0]
    rd = assnt.shape[1]
    pop_eq = propagate(TPM, int(t_max), dim)
    TPM_d1, TPM = get_M(TPM, dt, ms_lt)
    PAD = get_projector(assnt, pop_eq)
    projector = [assnt, PAD]
    #omega = np.dot(np.dot(assnt.T, TPM_d1), PAD)

    # calculate K1(t), K3(t), and C^EXACT(t)
    exact, ITS = dyn_with_ITS(TPM, projector, tsteps, dt, timearray)

    return exact, ITS


def dyn_with_ITS1(TPM, projector, tsteps, dt, timearray):
    """
    This function only computes the dimensionality reduced dynamics
    obtained via direct DR.
    """
    assnt = projector[0]
    PAD = projector[1]
    dim = assnt.shape[0]
    rd = assnt.shape[1]

    dynamics = np.zeros([rd, rd, tsteps], float)
    Tn_loop = np.identity(dim, float)
    ITS = np.zeros([1, int(rd - 1), tsteps], float)
    #dynamics[:, :, 0] = np.dot(assnt.T, PAD)

    # get exact dynamics, memory kernel, and generator
    for k in range(0, tsteps):
        # propagating the many-state model
        Tn_loop = np.dot(Tn_loop, TPM)
        TPAD = np.dot(Tn_loop, PAD)
        dynamics[:, :, k] = np.dot(assnt.T, TPAD)
        time = timearray[k]
        ITS[:, :, k] = implied_time_scale(dynamics[:, :, k], time, rd)

    return dynamics, ITS


def dyn_with_ITS(TPM, projector, tsteps, dt, timearray):
    """
    This function only computes the dimensionality reduced dynamics
    obtained via direct DR.
    """
    assnt = projector[0]
    PAD = projector[1]
    dim = assnt.shape[0]
    rd = assnt.shape[1]

    dynamics = np.zeros([rd, rd, tsteps], float)
    Tn_loop = np.identity(dim, float)
    ITS = np.zeros([1, int(rd - 1), tsteps], float)
    dynamics[:, :, 0] = np.dot(assnt.T, PAD)

    # get exact dynamics, memory kernel, and generator
    for k in range(1, tsteps):
        # propagating the many-state model
        Tn_loop = np.dot(Tn_loop, TPM)
        TPAD = np.dot(Tn_loop, PAD)
        dynamics[:, :, k] = np.dot(assnt.T, TPAD)
        time = timearray[k]
        ITS[:, :, k] = implied_time_scale(dynamics[:, :, k], time, rd)

    return dynamics, ITS


def implied_time_scale(tau, time, dim):
    """
    This function returns ITS(t) for each of N = (dimension - 1) eigenvalues

    Inputs:
    1. tau     - the TPM in the reduced space (matrix)
    2. time    - the time-step of this particular propagation (scalar)
    3. dim - the dimension of the reduced space (scalar)

    Outputs:
    1. ITS - 1 x num array containing the implied time scale value at this time
             step for the N-1 eigenvalues (1 x num dimensional array)
    """
    # the number of eigenvalues not equal to 1
    num = dim - 1

    # get the eigenvalues of tau and sort them
    eig, U = np.linalg.eig(tau)
    eig = np.sort(np.abs(eig), axis=None)
    #eig = np.sort(eig, axis=None)

    # loop to compute and store the ITS values
    ITS = np.zeros(num, float)
    for j in range(0, num):
        ITS[j] = - time / np.log(eig[j])

    return ITS


def psuedo_ITS(U, dim):
    """
    This function returns ITS(t) for each of N = (dimension - 1) eigenvalues

    Inputs:
    1. tau     - the TPM in the reduced space (matrix)
    2. time    - the time-step of this particular propagation (scalar)
    3. dim - the dimension of the reduced space (scalar)

    Outputs:
    1. ITS - 1 x num array containing the implied time scale value at this time
             step for the N-1 eigenvalues (1 x num dimensional array)
    """
    # the number of eigenvalues not equal to 1
    num = dim - 1

    # get the eigenvalues of tau and sort them
    eig, U = np.linalg.eig(U)
    eig = np.sort(np.abs(eig), axis=None)
    #eig = np.sort(eig, axis=None)

    # loop to compute and store the ITS values
    ITS = np.zeros(num, float)
    for j in range(0, num):
        ITS[j] = eig[j]

    return ITS


def diagonalize(matrix, dt):
    """
    This function diagonalizes a matrix and produces a propagator: e^{-Adt}
    """
    eig, P = np.linalg.eig(matrix)
    P_inv = np.linalg.inv(P)
    exp = np.exp(eig * dt)
    diag = np.diag(exp)
    rot = np.dot(P, np.dot(diag, P_inv))
    return rot


def ITS(TCM_file, assnt, num_lines, sym=True):
    dim, rd = assnt.shape
    time_array = [0.0]
    data = np.zeros([rd, rd, num_lines], float)
    data[:, :, 0] = np.identity(rd, float)

    pi_matrix = np.zeros([dim, dim, num_lines], float)
    projectors = []
    pop_eqs = []

    ms_ITS = np.zeros([1, int(dim - 1), num_lines], float)
    rd_ITS = np.zeros([1, int(rd - 1), num_lines], float)

    # get the PAD component of the projector
    C_temp = get_count_matrix(TCM_file, num_lines-1)
    C_temp = (C_temp + C_temp.T) / 2.0
    TPM_temp = normalize_columns(C_temp)
    pop_eq = propagate(TPM_temp, num_lines, dim)
    PAD = get_projector(assnt, pop_eq)


    for k in range(1, num_lines):
        time_array.append(float(k))

        # obtain count matrix and construct projectors and derivates
        C = get_count_matrix(tcm_file, int(ms_lt)-1)

        # symmetrize the count matrix
        if sym is True:
            C = (C + C.T) / 2.0

        TPM = normalize_columns(C)
        ms_ITS[:, :, k] = implied_time_scale(TPM, k, dim)

        # get data
        data[:, :, k] = np.dot(assnt.T, np.dot(TPM, PAD))
        pi_matrix[:, :, k] = np.diag(pop_eq)
        projectors.append(PAD)
        pop_eqs.append(pi_matrix)
        rd_ITS[:, :, k] = implied_time_scale(data[:, :, k], k, rd)
    time_array = np.array(time_array, float)
    
    return rd_ITS, ms_ITS


def gen_mem(tsteps, dt, dim_red, auxk1, auxk3):
    """
    tsteps - number of steps
    dt - time-step
    dim_red - reduced space dimension
    auxk1 - auxiliary kernel
    auxk3 - auxiliary kernel
    """
    mem = np.zeros_like(auxk1, float)
    temp = np.identity(dim_red, float) - dt * auxk3[:, :, 0] / 2.0
    prefac = np.linalg.inv(temp)

    mem_zero = auxk1[:, :, 0]
    mem[:, :, 0] = mem_zero

    for k in range(1, tsteps):
        value = np.zeros([dim_red, dim_red], float)

        for j in range(1, k):
            value += np.dot(auxk3[:, :, k - j], mem[:, :, j])
        value *= dt
        value += auxk1[:, :, k] + np.dot(auxk3[:, :, k], mem_zero) * dt / 2.0
        mem[:, :, k] = np.dot(prefac, value)
    
    return mem


def gme_evolution(mem, tau_0, omega, tsteps, cutoff, alpha, dt):

    #gme = np.zeros_like(mem, float)
    rd = omega.shape[0]
    gme = np.zeros([rd, rd, tsteps], float)
    # define some variables
    hdt = dt * 0.5
    #tc_max = int(cutoff / dt) + 1
    tc_max = cutoff

    # make the trapezoidal rule more convenient
    gme[:, :, 0] = tau_0
    prefac = omega + alpha * hdt * mem[:, :, 0]

    # get tau(1)
    k1 = np.dot(omega, tau_0)
    gme_tilde = tau_0 + dt * k1
    temp = hdt * np.dot(mem[:, :, 1], tau_0)
    k2 = np.dot(prefac, gme_tilde) + alpha * temp
    gme[:, :, 1] = tau_0 + hdt * (k1 + k2)

    # get tau(2)
    k1 = np.dot(prefac, gme[:, :, 1]) + alpha * temp
    gme_tilde = gme[:, :, 1] + dt * k1
    temp = dt * (np.dot(mem[:, :, 1], gme[:, :, 1]) +
                 np.dot(mem[:, :, 2], tau_0) / 2.0)
    k2 = np.dot(prefac, gme_tilde) + alpha * temp
    gme[:, :, 2] = gme[:, :, 1] + hdt * (k1 + k2)

    # generating tau(n+1) = tau(n) + hdt * (k1 + k2)
    for i in range(2, tsteps-1):

        # recycle ctemp from the previous step
        k1 = np.dot(prefac, gme[:, :, i]) + alpha * temp
        gme_tilde = gme[:, :, i] + dt * k1

        temp = np.zeros_like(tau_0)

        # set upper bound for the cutoff times
        if (i < tc_max):
            ub = i + 1
            temp += np.dot(mem[:, :, i + 1], tau_0) / 2.0

        else:
            ub = tc_max

        # trapezoidal rule
        for k in range(1, ub):
            temp += np.dot(mem[:, :, k], gme[:, :, i - k + 1])

        # make sure every term gets a dt
        temp *= dt

        k2 = np.dot(prefac, gme_tilde) + alpha * temp
        gme[:, :, i + 1] = gme[:, :, i] + hdt * (k1 + k2)

    return gme


def TTM(T, ub):
    dim = T.shape[0]
    C = np.zeros_like(T)

    # get the first two
    C[:, :, 0] = np.identity(dim)
    C[:, :, 1] = T[:, :, 1]

    # get the rest of them
    for k in range(2, T.shape[-1]):
        temp = np.zeros([dim, dim])

        for j in range(1, min(k, ub)):
            temp += np.dot(T[:, :, j], C[:, :, k-j])

        C[:, :, k] = temp + T[:, :, k]

    return C


def get_T(C):
    T = np.zeros_like(C)

    # Get the first two TTs
    T[:, :, 1] = C[:, :, 1]
    T[:, :, 2] = C[:, :, 2] - np.dot(T[:, :, 1], C[:, :, 1])

    # get the rest of em
    for k in range(3, C.shape[-1]):
        temp = np.zeros([C.shape[0], C.shape[0]])

        for j in range(1, k):
            temp += np.dot(T[:, :, j], C[:, :, k-j])

        T[:, :, k] = C[:, :, k] - temp

    return T

