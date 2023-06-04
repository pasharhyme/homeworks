import numpy as np
import matplotlib.pyplot as plt
import PVTProps as pvt

#Constants
SN  = [2, 4, 4, 2, 0, 4, 4] #Student number
GG  = 0.65 + (SN[6] + 1) / 50 #Total gas gravity
CO2 = (SN[4] + 1) / 200 #CO2 mole fraction
H2S = (SN[5] + 1) / 200 #H2S mole fraction
N2  = (SN[3] + 1) / 200 #N2 mole fraction
T   = np.array([60, 200, 400]) #F Temperature
PR  = np.array([15, 10000]) #psi Pressure range
P   = np.linspace(PR[0], PR[1], 100) #psi Pressure vector
TR  = np.array([(T + 460),] * np.size(P)).T #R Temperature

API  = 20 + SN[6]*3 #API Gravity
RSB  = (API*6)**1.2 + (SN[5]-5)*20 #Rsb
PGG  = 0.65 + (SN[6]+1)/50 #Produced gas gravity
TO   = 120 + SN[3]*3 #Oil temperature
TSEP = 120 #Seperator temperature
PSEP = 125 #Seperator pressure

CH4 = np.array([((SN[5]+1)*7/100), (-116.67+460), 666.4, 0.0104]) #CH4 properties
C2H6 = np.array([((SN[6]+1)*2/100), (89.92+460), 705.5, 0.0979]) #C2H4 properties
C3H8 = np.array([(1.0 - CH4[0] - C2H6[0]), (205.05 + 460), 616.0, 0.1522]) #C3H8 properties
R = 10.73 #Universal gas constant

#Correction methods
sCorrCriticProp = "Standing"
sCorrNonHC      = "WichertAziz"
sCorrViscosity  = "LeeGonzalesEakin"
sCoorelation    = "VasquezBeggs"

#Calculated pseudo critical properties
tpc = pvt.Tpc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
ppc = pvt.Ppc(GG, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)

#Calculated pseudo reduced critical properties 
trpc = TR / tpc
prpc = P / ppc

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

def zfactor_srk_eos():
    sT = 200 + 460 #Temperature R
    sP = np.array([15, 1000, 2000, 3000, 4000, 5000, 6000]) #Pressure psi
    stpc = (CH4[0]*CH4[1]) + (C2H6[0]*C2H6[1]) + (C3H8[0]*C3H8[1])
    sppc = (CH4[0]*CH4[2]) + (C2H6[0]*C2H6[2]) + (C3H8[0]*C3H8[2])
    saf = (CH4[0]*CH4[3]) + (C2H6[0]*C2H6[3]) + (C3H8[0]*C3H8[3])
    sTr = sT / stpc
    
    sa = 0.42727 * R**2 * stpc**2 / sppc
    sb = 0.08664 * R * stpc / sppc
    salpha = (1 + (0.480 + (1.574*saf) - (0.176*saf**2)) * (1 - stpc**0.5))**2
    sA = salpha * sa * sP / (R**2 * sT**2)
    sB = sb * sP / (R * sT)
    
    c = np.zeros((4, np.size(sP)))
    c[0,:] = 1
    c[1,:] = -1
    c[2,:] = (sA - sB - sB**2)
    c[3,:] = -sA
    
    roots = np.zeros((3, np.size(sP)))
    for i in range(3):
        roots[:,i] = np.real(np.roots(c[:,i]))
    print(roots)
    
#Figure generation and cleanup of unnecessary subplots
FIGSIZE = 5 #Size of the figure    
fig, ax = plt.subplots(3, 4, figsize=(4*FIGSIZE, 3*FIGSIZE))
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