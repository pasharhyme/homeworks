#File to hold the information regarding to the components that this program can process
import numpy as np #Dependency => pip install numpy
import sys

#Component = ["Pc" - psi, "Tc" - R, "Vc" - ft3/lbmole, "OMEGA" - 1]
C = np.array([["methane", 666.4, -116.67+460.0, 0.0958*16.043, 0.010],
              ["ethane", 706.5, 89.92+460.0, 0.0783*30.070, 0.098],
              ["propane", 616.0, 206.06+460.0, 0.0727*44.097, 0.153],
              ["n-Butane", 550.6, 305.62+460, 0.0703*58.123, 0.199],
              ["n-Pentane", 488.5, 385.8+460, 0.0875*72.150, 0.251]],
              dtype=object)

def constants(c, t):
    while True:
        try:
            row = np.where(t[:,0] == c)[0]
        except ValueError:
            print("Invalid entry. Terminating program...")
            sys.exit()
        if len(row) == 0:
            print("Invalid entry. Terminating program...")
            sys.exit()
        else:
            break

    return t[row.tolist()[0], 1:]

A = np.array([["1", "srk_eos"],
              ["2", "vdw_eos"],
              ["3", "pr_eos"]],
              dtype=object)

def eos(x):
    global model

    while True:
        try:
            model = np.where(A[:,0] == x)[0]
        except IndexError:
            print("Invalid entry. Terminating program...")
            sys.exit()
        else:
            break

    return A[model.tolist()[0], 1]