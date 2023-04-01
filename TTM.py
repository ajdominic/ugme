import numpy as np
import msm_gme as m


def TTM(T, ub, tsteps=1000):
    dim = T.shape[0]
    #C = np.zeros_like(T, dtype=float)
    C = np.zeros([dim, dim, tsteps])

    # instead of enforcing cut-off, set other elems to 0
    T2 = np.zeros_like(C)
    T2[:, :, :T.shape[-1]] = np.copy(T)

    # get the first two
    C[:, :, 0] = np.identity(dim, dtype=float)
    C[:, :, 1] = T2[:, :, 1]

    # get the rest of them
    #for k in range(2, T.shape[-1]):
    for k in range(2, tsteps):
        temp = np.zeros([dim, dim], dtype=float)

        #for j in range(1, min(k, ub)):
        for j in range(1, k):
            temp += np.dot(T2[:, :, j], C[:, :, k-j])

        C[:, :, k] = temp + np.dot(T2[:, :, k], C[:, :, 0])

    return C


def TTMv2(T, tsteps, tk):
    """
    Inputs
    1. T - the TTM 
    2. tsteps - total number of steps
    3. tk - TTM cut-off index

    Outputs
    1. C -  dynamics as a function of time
    """
    dim = T.shape[0]
    C = np.zeros([dim, dim, tsteps])

    # get the first two
    C[:, :, 0] = np.identity(dim, dtype=float)
    C[:, :, 1] = T[:, :, 1]

    # get the rest of them
    #for k in range(2, T.shape[-1]):
    for k in range(2, tsteps):
        temp = np.zeros([dim, dim], dtype=float)

        #for j in range(1, min(k, ub)):
        for j in range(1, min(k, tk)):
            temp += np.dot(T[:, :, j], C[:, :, k-j])

        if k < tk:
            C[:, :, k] = temp + np.dot(T[:, :, k], C[:, :, 0])
        else:
            C[:, :, k] = temp

    return C


def get_T(C):
    T = np.zeros_like(C, dtype=float)
    C0_inv = np.linalg.inv(C[:, :, 0])

    # Get the first two TTs
    T[:, :, 1] = np.dot(C[:, :, 1], C0_inv)
    term1 = C[:, :, 2] - np.dot(T[:, :, 1], C[:, :, 1])
    T[:, :, 2] = np.dot(term1, C0_inv)

    # get the rest of em
    for k in range(3, C.shape[-1]):
        temp = np.zeros([C.shape[0], C.shape[0]], dtype=float)

        for j in range(1, k):
            temp += np.dot(T[:, :, j], C[:, :, k-j])

        T[:, :, k] = np.dot(C[:, :, k] - temp, C0_inv)

    return T


def normalize_C(C):
    C_new = np.zeros_like(C)
    C0_inv = np.linalg.inv(C[:, :, 0])
    C_new[:, :, 0] = np.identity(C.shape[0], float)
    for k in range(1, C.shape[-1]):
        C_new[:, :, k] = np.dot(C[:, :, k], C0_inv)
    return C_new


def shifted_C(exact, lb, h=2000):
    """
    Inputs:
    exact - exact dynamics
    lb - the lowerbound
    h - length of the dynamics

    Outputs:
    C - the dynamics
    T - the transfer matrix
    """
    ub = lb + h
    shifted_dynamics = exact[:, :, lb:ub]

    # get transfer matrix
    TTM_shifted = get_T(shifted_dynamics)

    # get dynamics from the transfer matrix
    #C_shifted = t.TTM(TTM_shifted, h)
    return TTM_shifted
    
    #return C_shifted, TTM_shifted


def tests(EXACT, val=10):
    """
    This function takes in the exact dynamics and computes
    two generic tests:
    test 1: does the TTM reproduce the exact dynamics?
    test 2: does the TTM truncation work?

    val sets the specific kernel cutoff time used
    Ex: when val = 20 for A2P, this corresponds to tk=5ps
    """
    tsteps = EXACT.shape[-1]
    tk = tsteps
    T_test1 = get_T(EXACT)
    C_test1 = TTMv2(T_test1, tsteps, tk)

    tk = int(tsteps/val)
    T_test2 = get_T(EXACT)
    C_test2 = TTMv2(T_test2, tsteps, tk)

    print(f"The error for test 1 is {m.rmse(EXACT, C_test1)}")
    print(f"The error for test 2 is {m.rmse(EXACT, C_test2)}")
