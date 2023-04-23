#Construct isotherms of a pure hydrocarbon component using Soave-Redlich-Kwong equation of state
import numpy as np
import matplotlib.pyplot as plt

#Constants for the component that will be studied
TC = 206.06 + 460 #R
PC = 616.0 #psi
VC = 44.097 * 0.0727 #ft3/lbmole
OMEGA = 0.153 #1
#Other Constants
V = 100 #Magnitude of the volume vector
NT = 6 #Number of isotherms to be worked with
T = np.array([60, 100, 140, 180, 206.6, 220]) #Temperature in F
NPS = 4 #Number of PSat Estimation values
PSE = np.array([108.83, 191.35, 311.86, 478.87]) #Estimated Saturaion Pressure values

R = 10.732 #psi*ft3/(lbmol*R)

tt = T + 460 #Temperature in R
tr = tt/TC #Reduced Temperature
v = np.linspace(1, V, V * 10) #Vector of Volume

#Soave-Redlich-Kwong Equation of State
def srk_eos(tx, trx, vx):
    p = np.zeros((NT, V * 10))
    a = 0.42727 * R**2 * TC**2 / PC
    b = 0.08664 * R * TC / PC
    m = 0.480 + 1.574 * OMEGA - 0.176 * OMEGA**2
    alpha = (1 + m * (1 - np.sqrt(trx)))**2
    for i in range(NT): #Calculate pressure for each isotherm
        p[i, :] = (R * tx[i] / (vx - b)) - (a * alpha[i] / (vx * (vx + b)))
    return p

#Cubic form of Soave-Redlich-Kwong Equation of State
def srk_eos_psat(ttx, trx):
    c = np.zeros((NPS, NPS)) #Array to store coefficients for every PSE
    a = 0.42727 * R**2 * TC**2 / PC
    b = 0.08664 * R * TC / PC
    m = 0.480 + 1.574 * OMEGA - 0.176 * OMEGA**2
    alpha = (1 + m * (1 - np.sqrt(trx)))**2
    for i in range(NPS): #Generate coefficients for each PSE
        c[i, 0] = PSE[i]
        c[i, 1] = -(R * ttx[i])
        c[i, 2] = (a * alpha[i]) - (b * R * ttx[i]) - (b**2 * PSE[i])
        c[i, 3] = -(a * b * alpha[i])
    roots = np.zeros((NPS, 3))
    for i in range(NPS): #Calculate roots of the equation for each PSE
        roots[i, :] = np.real(np.roots(c[i, :].T))
    roots = np.fliplr(roots)
    ai = np.zeros((NPS, 3)) #Roots of the integral of srk eos
    for i in range(NPS): #Generate I values
        ai[i, 0] = (R*ttx[i]*np.log(abs(roots[i, 0]-b))) + ((a*alpha[i]/b)*np.log(abs((b/roots[i,0])+1)))
        ai[i, 1] = (R*ttx[i]*np.log(abs(roots[i, 1]-b))) + ((a*alpha[i]/b)*np.log(abs((b/roots[i,1])+1)))
        ai[i, 2] = (R*ttx[i]*np.log(abs(roots[i, 2]-b))) + ((a*alpha[i]/b)*np.log(abs((b/roots[i,2])+1)))
    aj = np.zeros((NPS, 2)) #Area under Psat line
    for i in range(NPS):
        aj[i, 0] = (roots[i, 1] - roots[i, 0]) * PSE[i]
        aj[i, 1] = (roots[i, 2] - roots[i, 1]) * PSE[i]
    aa = np.zeros((NPS, 2)) #Area between loops and Psat line
    for i in range(NPS):
        aa[i, 0] = aj[i, 0] - (ai[i, 1] - ai[i, 0])
        aa[i, 1] = (ai[i, 2] - ai[i, 1]) - aj[i, 1]
    global adiff
    adiff = np.zeros((NPS, 1)) #Area difference should be 0 to obtain correct Psat values
    for i in range(NPS):
        adiff[i, 0] = aa[i, 0] - aa[i, 1]
    global PSEPlot
    PSEPlot = np.array([PSE,] * 3).T #Reformat estimation values to work with 3 roots
    return roots

def psat_corr(trx):
    trcorr = np.linspace(0.4, 1.0, 100 * 10)
    tcorr = (trcorr * TC) - 460
    a = (5.92714) - (6.09648/trcorr) - (1.28866*np.log(trcorr)) + (0.16934*(trcorr)**6)
    b = (15.2518) - (15.6875/trcorr) - (13.4721*np.log(trcorr)) + (0.4357*(trcorr)**6)
    Psatval = (np.exp(a + (OMEGA * b))) * PC
    Psat = np.vstack((tcorr, Psatval)).T
    return Psat

#Create a figure with 3 subplots and resize it for better visibility
fig, axes = plt.subplots(1, 3)
FSIZE = 5
fig.set_size_inches(FSIZE * 3, FSIZE)

#[Subplot 1] for pressure-volume graph
for i in range(NT):
    axes[0].semilogx(v, srk_eos(tt, tr, v)[i, :], label = f"{T[i]} °F")

#[Subplot 2] for pressure-volume graph with vapor pressure values
for i in range(NT):
    axes[1].semilogx(v, srk_eos(tt, tr, v)[i, :], label = f"{T[i]} °F")
for i in range(NPS):
    axes[1].semilogx(srk_eos_psat(tt, tr)[i, :], PSEPlot[i, :], "o:k")

#[Subplot 3] for comparing estimated and calculated vapor pressure values
axes[2].plot(psat_corr(tr)[:,0], psat_corr(tr)[:,1], label="PSat Correlated")
axes[2].scatter(T[0:np.size(PSE)], PSE, color="red", label="PSat Estimated")

#Set rules for first two subplots
for i in range(2):
    axes[i].set_box_aspect(1)
    axes[i].set_xlabel('Volume, ft^3')
    axes[i].set_ylabel('Pressure, psi')
    axes[i].set_xlim([1.1, V])
    axes[i].set_ylim([0, 800])
    axes[i].legend(loc = 1)

#Set rules for the third subplot
axes[2].set_box_aspect(1)
axes[2].set_xlabel('Temperature, °F')
axes[2].set_ylabel('Vapor Pressure, psi')
axes[2].legend(loc = 0)

plt.show() #Render the figure

#Report area difference as error
print("Estimation Error:", abs(adiff.T))