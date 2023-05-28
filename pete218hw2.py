import numpy as np
import matplotlib.pyplot as plt
import PVTProps as pvt

#Constants
SN  = [2, 4, 4, 2, 0, 4, 4] #Student number
GG  = 0.65 + (SN[6] + 1) / 50 #Total gas gravity
CO2 = (SN[4] + 1) / 200 #CO2 mole fraction
H2S = (SN[5] + 1) / 200 #H2S mole fraction
N2  = (SN[3] + 1) / 200 #N2 mole fraction
T  = np.array([60, 200, 400]) #F Temperature
PR  = np.array([15, 10000]) #psi Pressure range
P   = np.linspace(PR[0], PR[1], 1000) #psi Pressure vector (1 x 1000)
TR  = np.array([(T + 460),] * np.size(P)).T #R Temperature (3 x 1000)

API = 20 + SN[6]*3 #API Gravity
RSB = (API*6)**1.2 + (SN[5]-5)*20 #Rsb
PGG = 0.65 + (SN[6]+1)/50 #Produced gas gravity
TO = 120 + SN[3]*3 #Oil temperature
TSEP = 120 #Seperator temperature
PSEP = 125 #Seperator pressure

#Correction methods
sCorrCriticProp     = "Standing"
sCorrNonHC          = "WichertAziz"
sCorrViscosity      = "LeeGonzalesEakin"
sCoorelation        = "VasquezBeggs"

#Calculated pseudo critical properties
tpc = pvt.Tpc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
ppc = pvt.Ppc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)

#Calculated pseudo reduced critical properties 
trpc = TR / tpc #(3 x 1000)
prpc = P / ppc  #(1 x 1000)

#Z-Factor calculation
zFactor = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        zFactor[i,j] = pvt.ZFactor(trpc[i,j], prpc[j])

#Gas compressibility calculation
cg = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        cg[i,j] = pvt.Cg(T[i], P[j], GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)

#Gas viscosity calculation
ug = np.zeros((np.size(T), np.size(P)))
for i in range(np.size(T)):
    for j in range(np.size(P)):
        ug[i,j] = pvt.Ug(T[i], P[j], GG, N2, CO2, H2S, None, None, sCorrViscosity, sCorrCriticProp, sCorrNonHC)

#Solution Gas Oil Ratio calculation
rss = np.zeros((np.size(P)))
rsv = np.zeros((np.size(P)))
for i in range(np.size(P)):
    rss[i] = pvt.Rs(API, RSB, PGG, TO, P[i], "Standing", TSEP, PSEP)
    rsv[i] = pvt.Rs(API, RSB, PGG, TO, P[i], sCoorelation, TSEP, PSEP)

#Oil Formation Volume Factor calculation
bos = np.zeros((np.size(P)))
bov = np.zeros((np.size(P)))
for i in range(np.size(P)):
    bos[i] = pvt.Bo(API, RSB, PGG, TO, P[i], "Standing", TSEP, PSEP)
    bov[i] = pvt.Bo(API, RSB, PGG, TO, P[i], sCoorelation, TSEP, PSEP)

#Oil Viscosity calculation
uo = np.zeros((np.size(P)))
for i in range(np.size(P)):
    uo[i] = pvt.Uo(API, RSB, PGG, TO, P[i])

#Oil Compressibility calculation
cos = np.zeros((np.size(P)))
cov = np.zeros((np.size(P)))
for i in range(np.size(P)):
    cos[i] = pvt.Co(API, RSB, PGG, TO, P[i], "Standing", TSEP, PSEP)
    cov[i] = pvt.Co(API, RSB, PGG, TO, P[i], sCoorelation, TSEP, PSEP)
 
#Figure generation and cleanup of unnecessary subplots    
fig, ax = plt.subplots(3, 4, figsize=(4*10, 3*10))
fig.delaxes(ax[0][3])
fig.delaxes(ax[2][1])
fig.delaxes(ax[2][2])
fig.delaxes(ax[2][3])

#Plots for Question 1
for i in range(np.size(T)):
    ax[0,0].semilogx(P, zFactor[i,:], label=f"{T[i]} °F")
    ax[0,0].legend()
for i in range(np.size(T)):
    ax[0,1].semilogx(P, cg[i,:], label=f"{T[i]} °F")
    ax[0,1].legend()
for i in range(np.size(T)):
    ax[0,2].semilogx(P, ug[i,:], label=f"{T[i]} °F")
    ax[0,2].legend()

#Plots for Question 2
ax[1,0].semilogx(P, rss, label="Standing")
ax[1,0].semilogx(P, rsv, label="Vasquez-Beggs")
ax[1,0].legend()

ax[1,1].semilogx(P, bos, label="Standing")
ax[1,1].semilogx(P, bov, label="Vasquez-Beggs")
ax[1,1].legend()

ax[1,2].semilogx(P, uo)

ax[1,3].semilogx(P, cos, label="Standing")
ax[1,3].semilogx(P, cov, label="Vasquez-Beggs")
ax[1,3].legend()

plt.show() #Render the figure