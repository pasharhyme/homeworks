import numpy as np

##############Gas Props##############

"""
T = temperature of interest (degrees F)
P = pressure of interest (psia)
Gg = specific gravity of gas including non-hydrocarbon components
N2 = % Nitrogen
CO2 = % CO2
H2S = % H2S
"""
#Calculation of Critical Temperature for Miscelaneous Gases
def Tpc(Gg, N2, CO2, H2S, sCorrCriticProp="Standing", sCorrNonHC="WichertAziz"):
    yN2 = N2 / 100
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    cGg = (Gg - 0.9672 * yN2 - 1.5197 * yCO2 - 1.1768 * yH2S) / (1 - yN2 - yCO2 - yH2S)
    if (sCorrCriticProp=="Standing"):
        Tpc = 168 + 325 * cGg - 12.5 * cGg ** 2
    elif sCorrCriticProp=="Sutton":
        Tpc = 169.2 + 349.5 * cGg - 74 * cGg ** 2
    else:
        raise ValueError("Invalid sCorrCriticProp value")
        return None
    Tpc = (1 - yN2 - yCO2 - yH2S) * Tpc + 227.5 * yN2 + 547.9 * yCO2 + 672.4 * yH2S
    if (sCorrNonHC=="WichertAziz"):
        Cwa = 120 * ((yCO2 + yH2S) ** 0.9 - (yCO2 + yH2S) ** 1.6) + 15 * (yH2S ** 0.5 - yH2S ** 4)
        Tpc = Tpc - Cwa
    elif (sCorrNonHC=="CarrKobayashiBurrows"):
        Tpc = Tpc - 80 * yCO2 + 130 * yH2S - 250 * yN2
    else:
        raise ValueError("Invalid sCorrNonHC value")
        return None
    return Tpc

#Correct critical temperature
def CorrTpc(Tc, N2, CO2, H2S, sCorrNonHC="WichertAziz"):
    yN2 = N2 / 100
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    if (sCorrNonHC=="WichertAziz"):
        Cwa = 120 * ((yCO2 + yH2S) ** 0.9 - (yCO2 + yH2S) ** 1.6) + 15 * (yH2S ** 0.5 - yH2S ** 4)
        CorrTpc = Tc - Cwa
    elif (sCorrNonHC=="CarrKobayashiBurrows"):
        CorrTpc = Tc - 80 * yCO2 + 130 * yH2S - 250 * yN2
    else:
        raise ValueError("Invalid sCorrNonHC value")
        return None
    return CorrTpc

#Correction coefficient for Wichert-Aziz method
def WichAzizCorrTerm(CO2, H2S):
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    WichAzizCorrTerm = 120 * ((yCO2 + yH2S) ** 0.9 - (yCO2 + yH2S) ** 1.6) + 15 * (yH2S ** 0.5 - yH2S ** 4)
    return WichAzizCorrTerm

#Calculation of Critical Pressure for Miscelaneous Gases
def Ppc(Gg, N2, CO2, H2S, sCorrCriticProp="Standing", sCorrNonHC="WichertAziz"):
    yN2 = N2 / 100
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    cGg = (Gg - 0.9672 * yN2 - 1.5197 * yCO2 - 1.1768 * yH2S) / (1 - yN2 - yCO2 - yH2S)
    if (sCorrCriticProp=="Standing"):
        Ppc = 677 + 15 * cGg - 37.5 * cGg ** 2
    elif (sCorrCriticProp=="Sutton"):
        Ppc = 756.8 - 131 * cGg - 3.6 * cGg ** 2
    else:
        raise ValueError("Invalid sCorrCriticProp value")
        return None
    Ppc = (1 - yN2 - yCO2 - yH2S) * Ppc + 493.1 * yN2 + 1071 * yCO2 + 1306 * yH2S
    if (sCorrNonHC=="WichertAziz"):
        Tc = 168 + 325 * cGg - 12.5 * cGg ** 2
        Tc = (1 - yN2 - yCO2 - yH2S) * Tc + 227.5 * yN2 + 547.9 * yCO2 + 672.4 * yH2S
        Cwa = 120 * ((yCO2 + yH2S) ** 0.9 - (yCO2 + yH2S) ** 1.6) + 15 * (yH2S ** 0.5 - yH2S ** 4)
        cTc = Tc - Cwa
        Ppc = Ppc * cTc / (Tc + yH2S * (1 - yH2S) * Cwa)
    elif (sCorrNonHC=="CarrKobayashiBurrows"):
        Ppc = Ppc + 440 * yCO2 + 600 * yH2S - 170 * yN2
    else:
        raise ValueError("Invalid sCorrNonHC value")
        return None
    return Ppc

#Correct critical pressure
def CorrPpc(Pc, Tc, N2, CO2, H2S, sCorrNonHC="WichertAziz"):
    yN2 = N2 / 100
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    if (sCorrNonHC=="WichertAziz"):
        Cwa = 120 * ((yCO2 + yH2S) ** 0.9 - (yCO2 + yH2S) ** 1.6) + 15 * (yH2S ** 0.5 - yH2S ** 4)
        cTc = Tc - Cwa
        CorrPpc = Pc * cTc / (Tc + yH2S * (1 - yH2S) * Cwa)
    elif (sCorrNonHC=="CarrKobayashiBurrows"):
        CorrPpc = Pc + 440 * yCO2 + 600 * yH2S - 170 * yN2
    else:
        raise ValueError("Invalid sCorrNonHC value")
        return None
    return CorrPpc

#Calculation of Z Factor for Miscelaneous Gases using Pr and Tr
def ZFactor(Tr, Pr):
    a = 0.064225133
    b = 0.53530771 * Tr - 0.61232032
    c = 0.31506237 * Tr - 1.0467099 - 0.57832729 / Tr ** 2
    d = Tr
    e = 0.68157001 / Tr ** 2
    f = 0.68446549
    g = 0.27 * Pr
    """
    #Newton-Raphson iteration commented due to unstability and convergence problems
    rho = 0.27 * Pr / Tr  #Initial guess
    rhoold = rho
    for i in range(100):
        frho = a * rho ** 6 + b * rho ** 3 + c * rho ** 2 + d * rho + e * rho ** 3 * (1 + f * rho ** 2) * np.exp(-f * rho ** 2) - g
        dfrho = 6 * a * rho ** 5 + 3 * b * rho ** 2 + 2 * c * rho + d + e * rho ** 2 * (3 + f * rho ** 2 * (3 - 2 * f * rho ** 2)) * np.exp(-f * rho ** 2)
        rho = rho - frho / dfrho
        test = abs((rho - rhoold) / rho)
        if test < 0.00001:
            break
        rhoold = rho
    ZFactor = 0.27 * Pr / rho / Tr
    """

    #Find the solution by bi-section method
    #Initial interval of Z-Factor
    Za = 0.0001
    Zb = 4
    for i in range(100):
        Zx = (Za + Zb) / 2
        rho = 0.27 * Pr / (Zx * Tr)
        frho = a * rho ** 6 + b * rho ** 3 + c * rho ** 2 + d * rho + e * rho ** 3 * (1 + f * rho ** 2) * np.exp(-f * rho ** 2) - g
        if (frho > 0):
            Za = Zx
        else:
            Zb = Zx
        if (abs(Za - Zb) < 0.0001):
            break
    ZFactor = 0.27 * Pr / rho / Tr
    return ZFactor

#Calculation of Gas Formation Volume Factor for Miscelaneous Gases
def Bg(T, P, Gg, N2, CO2, H2S, sCorrCriticProp="Standing", sCorrNonHC="WichertAziz"):
    Tc = Tpc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Pc = Ppc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Tr = (T + 460) / Tc
    Pr = P / Pc
    Zm = ZFactor(Tr, Pr)
    Bg = 0.02827 * Zm * (T + 460) / P
    return Bg

#Calculation of pseudoreduced gas compressibility
def prCg(Tr, Pr):
    a = 0.064225133
    b = 0.53530771 * Tr - 0.61232032
    c = 0.31506237 * Tr - 1.0467099 - 0.57832729 / Tr ** 2
    d = Tr
    e = 0.68157001 / Tr ** 2
    f = 0.68446549
    g = 0.27 * Pr
    Zm = ZFactor(Tr, Pr)
    rho = 0.27 * Pr / Zm / Tr
    der = 1 / rho / Tr * (5 * a * rho ** 5 + 2 * b * rho ** 2 + c * rho + 2 * e * rho ** 2 * (1 + f * rho ** 2 - f ** 2 * rho ** 4) * np.exp(-f * rho ** 2))
    prCg = 1 / Pr / (1 + rho / Zm * der)
    return prCg

#Calculation of Gas Compressibility for Miscelaneous Gases
def Cg(T, P, Gg, N2, CO2, H2S, sCorrCriticProp="Standing",sCorrNonHC="WichertAziz"):
    Tc = Tpc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Pc = Ppc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Tr = (T + 460) / Tc
    Pr = P / Pc
    Cr = prCg(Tr, Pr)
    Cg = Cr / Pc
    return Cg

def Ug(T, P, Gg, N2, CO2, H2S, uTpc=None,uPpc=None,sCorrViscosity="LeeGonzalesEakin",sCorrCriticProp="Standing",sCorrNonHC="WichertAziz"):
    yN2 = N2 / 100
    yCO2 = CO2 / 100
    yH2S = H2S / 100
    if uTpc == None:
        Tc = Tpc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    else:
        Tc = uTpc
    if uPpc == None:
        Pc = Ppc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    else:
        Pc = uPpc
    Tr = (T + 460) / Tc
    Pr = P / Pc
    if (sCorrViscosity=="LeeGonzalesEakin"):
        Zm = ZFactor(Tr, Pr)
        MW = 28.97 * Gg
        a = (9.4 + 0.02 * MW) * (T + 460) ** 1.5 / (209 + 19 * MW + (T + 460)) / 10000
        b = 3.5 + 986 / (T + 460) + 0.01 * MW
        c = 2.4 - 0.2 * b
        rho = P * MW / Zm / 10.73 / (T + 460)
        Ug = a * np.exp(b * (rho / 62.4) ** c)
    elif (sCorrViscosity=="CarrKobayashiBurrows"):
        a0 = -2.4621182
        a1 = 2.97054714
        a2 = -0.286264054
        a3 = 0.008054205
        a4 = 2.80860949
        a5 = -3.49803305
        a6 = 0.36037302
        a7 = -0.01044324
        a8 = -0.793385684
        a9 = 1.39643306
        a10 = -0.149144925
        a11 = 0.004410155
        a12 = 0.083938718
        a13 = -0.186408848
        a14 = 0.020336788
        a15 = -0.000609579
        logGg = np.log(Gg) / np.log(10)
        uN2corr = (0.00959 + 0.00848 * logGg) * yN2
        uCO2corr = (0.00624 + 0.00908 * logGg) * yCO2
        uH2Scorr = (0.00373 + 0.00849 * logGg) * yH2S
        u1Atm = 0.008188 - 0.00615 * logGg + (0.00001709 - 0.000002062 * Gg) * T
        u1AtmCorr = u1Atm + uN2corr + uCO2corr + uH2Scorr
        u1UgP1 = a0 + a1 * Pr + a2 * Pr ** 2 + a3 * Pr ** 3
        u1UgP2 = (a4 + a5 * Pr + a6 * Pr ** 2 + a7 * Pr ** 3) * Tr
        u1UgP3 = (a8 + a9 * Pr + a10 * Pr ** 2 + a11 * Pr ** 3) * Tr ** 2
        u1UgP4 = (a12 + a13 * Pr + a14 * Pr ** 2 + a15 * Pr ** 3) * Tr ** 3
        u1UgS = u1UgP1 + u1UgP2 + u1UgP3 + u1UgP4
        u1Ug = np.exp(u1UgS) / Tr
        Ug = u1AtmCorr * u1Ug
    else:
        raise ValueError("Invalid sCorrViscosity value")
        return None
    return Ug

#Calculation of Gas Psuedopressure for Miscelaneous Gases
def mP(T, P, Gg, N2, CO2, H2S, sCorrViscosity="LeeGonzalesEakin",sCorrCriticProp="Standing",sCorrNonHC="WichertAziz"):
    mP = 0
    Pold = 0
    Xold = 0
    Pstep = P / 20
    Tc = Tpc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Pc = Ppc(Gg, N2, CO2, H2S, sCorrCriticProp, sCorrNonHC)
    Tr = (T + 460) / Tc
    for i in range(20):
        Pnew = Pold + Pstep
        Pr = Pnew / Pc
        Xnew = 2 * Pnew / ZFactor(Tr, Pr) / Ug(T, Pnew, Gg, N2, CO2, H2S, None, None, sCorrViscosity, sCorrCriticProp, sCorrNonHC)
        mP = mP + (Xold + Xnew) / 2 * Pstep
        Pold = Pnew
        Xold = Xnew
    return mP

##############Oil Props##############

"""
API = crude API gravity
Rsb = solution gas/oil ratio at bubble point pressure (SCF/STBO)
T = temperature of interest (degrees F)
Gg = gas specific gravity (air=1)
Tsep = separator temperature (degrees F)
Psep = separator pressure (psia)
"""
#Calculation of bubble point pressure
def BPP(API, Rsb, Gg, T, sCorrelation="Standing", Tsep=None, Psep=None):
    if (sCorrelation == "Standing"):
        BPP = 18.2 * ((Rsb / Gg) ** 0.83 * np.float64(10) ** (0.00091 * T - 0.0125 * API) - 1.4)
    elif (sCorrelation == "VasquezBeggs"):
        if Tsep==None or Psep==None:
            raise ValueError("Invalid sCorrelation value")
            return
        sGg = Gg * (1 + 0.00005912 * API * Tsep * 0.434294 * np.log(Psep / 114.7))
        if API > 30:
            C1 = 0.0178
            C2 = 1.187
            C3 = 23.931
        else:
            C1 = 0.0362
            C2 = 1.0937
            C3 = 25.724
        BPP = (Rsb / (C1 * sGg * np.exp(C3 * API / (460 + T)))) ** (1 / C2)
    else:
        raise ValueError("Invalid sCorrelation value")
        return None
    return BPP
#Calculation of solution gas/oil ratio at specified temperature and pressure for the given field data (SCF/STBO)
def Rs(API, Rsb, Gg, T, P, sCorrelation= "Standing", Tsep=None, Psep=None):
    Pb = BPP(API, Rsb, Gg, T, sCorrelation, Tsep, Psep)
    if P <= Pb:
        Rs = xRs(API, Gg, T, P, sCorrelation, Tsep, Psep)
    else:
        Rs = Rsb
    return Rs
#Calculation of solution gas/oil ratio at specified temperature and pressure for the given field data (SCF/STBO)
def xRs(API, Gg, T, P, sCorrelation="Standing", Tsep=None, Psep=None):
    if (sCorrelation == "Standing"):
        xRs = ((P / 18.2 + 1.4) / (np.float64(10) ** (0.00091 * T - 0.0125 * API))) ** (np.float64(1) / 0.83) * Gg
    elif (sCorrelation == "VasquezBeggs"):
        if Tsep==None or Psep==None:
            raise ValueError("Invalid sCorrelation value")
            return None
        sGg = Gg * (1 + 0.00005912 * API * Tsep * 0.434294 * np.log(Psep / 114.7))
        if API > 30:
            C1 = 0.0178
            C2 = 1.187
            C3 = 23.931
        else:
            C1 = 0.0362
            C2 = 1.0937
            C3 = 25.724
        xRs = C1 * sGg * P ** C2 * np.exp(C3 * (API / (T + 460)))
    else:
        raise ValueError("Invalid sCorrelation value")
        return None
    return xRs
#Calculation of oil formation volume factor
def Bo(API, Rsb, Gg, T, P, sCorrelation="Standing", Tsep=None, Psep=None):
    Pb = BPP(API, Rsb, Gg, T, sCorrelation, Tsep, Psep)
    if (P <= Pb):
        GOR = xRs(API, Gg, T, P, sCorrelation, Tsep, Psep)
    else:
        GOR = Rsb
    if (sCorrelation == "Standing"):
        Go = 141.5 / (API + 131.5)
        CBob = GOR * (Gg / Go) ** 0.5 + 1.25 * T
        Bo = 0.9759 + 0.00012 * CBob ** 1.2
        if P > Pb:     #correct if pressure is greater than bubble point
            xCo = Co(API, Rsb, Gg, T, P, sCorrelation, Tsep, Psep)
            Bo = Bo * np.exp(xCo * (Pb - P))
    elif (sCorrelation == "VasquezBeggs"):
        sGg = Gg * (1 + 0.00005912 * API * Tsep * 0.434294 * np.log(Psep / 114.7))
        if API > 30:
            C1 = 0.000467
            C2 = 0.000011
            C3 = 0.000000001337
        else:
            C1 = 0.0004677
            C2 = 0.00001751
            C3 = -0.00000001811
        Bo = 1 + C1 * GOR + C2 * (T - 60) * (API / sGg) + C3 * GOR * (T - 60) * (API / sGg)
        if P > Pb:     #correct if pressure is greater than bubble point
            xCo = Co(API, Rsb, Gg, T, P, sCorrelation, Tsep, Psep)
            Bo = Bo * np.exp(xCo * (Pb - P))
    else:
        raise ValueError("Invalid sCorrelation value")
        return None
    return Bo
#Calculation of oil viscosity
def Uo(API, Rsb, Gg, T, P):
    Pb = BPP(API, Rsb, Gg, T)
    if (P <= Pb):
        GOR = xRs(API, Gg, T, P)
    else:
        GOR = Rsb
    c = 3.0324 - 0.02023 * API
    b = 10 ** c
    a = b * T ** -1.163
    Uod = 10 ** a - 1
    a = 10.715 * (GOR + 100) ** -0.515
    b = 5.44 * (GOR + 150) ** -0.338
    Uobo = a * Uod ** b
    if P <= Pb:
        Uo = a * Uod ** b
    else:
        Uobp = a * Uod ** b
        a = 2.6 * P ** 1.187 * np.exp(-0.0000898 * P - 11.513)
        Uo = Uobp * (P / Pb) ** a
    return Uo
#Calculation of oil compressibility above bubble point only
def Co(API, Rsb, Gg, T, P, sCorrelation="Standing", Tsep=None, Psep=None):
    Pb = BPP(API, Rsb, Gg, T, sCorrelation, Tsep, Psep)
    if (P <= Pb):
        dP = P * 0.01
        P1 = P - dP / 2
        P2 = P + dP / 2
        if P2 > Pb:
            P2 = Pb
            P1 = Pb - dP
        Bo1 = Bo(API, Rsb, Gg, T, P1, sCorrelation, Tsep, Psep)
        Bo2 = Bo(API, Rsb, Gg, T, P2, sCorrelation, Tsep, Psep)
        dBo = (Bo2 - Bo1) / dP
        Rs1 = Rs(API, Rsb, Gg, T, P1, sCorrelation, Tsep, Psep)
        Rs2 = Rs(API, Rsb, Gg, T, P2, sCorrelation, Tsep, Psep)
        dRs = (Rs2 - Rs1) / dP
        BoP = (Bo1 + Bo2) / 2
        BgP = Bg(T, P, Gg, 0, 0, 0)
        Co = (-1 / BoP) * (dBo - BgP * dRs)
        return Co
    if (sCorrelation == "Standing"):
        sGg = Gg
    elif (sCorrelation == "VasquezBeggs"):
        sGg = Gg * (1 + 0.00005912 * API * Tsep * 0.434294 * np.log(Psep / 114.7))
    else:
        raise ValueError("Invalid sCorrelation value")
        return None
    a1 = -1433
    a2 = 5
    a3 = 17.2
    a4 = -1180
    a5 = 12.61
    a6 = 10 ** 5
    Co = (a1 + a2 * Rsb + a3 * T + a4 * sGg + a5 * API) / (a6 * P)
    return Co

##############Water Props##############

"""
Rsw = gas/water ratio (SCF/STBW)
T = temperature of interest (degrees F)
P = pressure (psia)
Csalt = Csalt concentration (weight %)
"""
def Cw(T, P, Rsw, Csalt):
    a = 3.8546 - 0.000134 * P
    b = -0.01052 + 0.000000477 * P
    c = 0.000039267 - 0.00000000088 * P
    Cw = (a + b * T + c * T ** 2) / np.float64(1000000)
    Cw = Cw * (1 + 0.0089 * Rsw)   #Dissolved gas correction
    Cw = Cw * ((-0.052 + 0.00027 * T - 0.00000114 * T ** 2 + 0.000000001121 * T ** 3) * Csalt ** 0.7 + 1)
    return Cw
def Bw(T, P, Csalt):  #Gas saturated water
    a = 0.9911 + 0.0000635 * T + 0.00000085 * T ** 2
    b = -0.000001093 - 0.000000003497 * T + 0.00000000000457 * T ** 2
    c = -0.00000000005 + 6.429E-13 * T - 1.43E-15 * T ** 2
    Bw = a + b * P + c * P ** 2
    Bw = Bw * ((0.000000051 * P + (0.00000547 - 0.000000000195 * P) * (T - 60) + (-0.0000000323 + 0.00000000000085 * P) * (T - 60) ** 2) * Csalt + 1)
    return Bw
def Uw(T, P, Csalt):
    Tc = 5 / 9 * (T - 32)
    tk = Tc + 273.15
    Sum = -7.419242 * (0.65 - 0.01 * Tc) ** (1 - 1)
    Sum = Sum - 0.29721 * (0.65 - 0.01 * Tc) ** (2 - 1)
    Sum = Sum - 0.1155286 * (0.65 - 0.01 * Tc) ** (3 - 1)
    Sum = Sum - 0.008685635 * (0.65 - 0.01 * Tc) ** (4 - 1)
    Sum = Sum + 0.001094098 * (0.65 - 0.01 * Tc) ** (5 - 1)
    Sum = Sum + 0.00439993 * (0.65 - 0.01 * Tc) ** (6 - 1)
    Sum = Sum + 0.002520658 * (0.65 - 0.01 * Tc) ** (7 - 1)
    Sum = Sum + 0.0005218684 * (0.65 - 0.01 * Tc) ** (8 - 1)
    psat = 22088 * np.exp((374.136 - Tc) * Sum / tk)
    Uw = 0.02414 * 10 ** (247.8 / (tk - 140)) * (1 + (P / 14.504 - psat) * 0.0000010467 * (tk - 305))
    Uw = Uw * (1 - 0.00187 * Csalt ** 0.5 + 0.000218 * Csalt ** 2.5 + (T ** 0.5 - 0.0135 * T) * (0.00276 * Csalt - 0.000344 * Csalt ** 1.5))
    return Uw
def RSwat(T, P, Csalt):
    a = 2.12 + 0.00345 * T - 0.0000359 * T ** 2
    b = 0.0107 - 0.0000526 * T + 0.000000148 * T ** 2
    c = -0.000000875 + 0.0000000039 * T - 0.0000000000102 * T ** 2
    RSwat = a + b * P + c * P ** 2
    RSwat = RSwat * (1 - (0.0753 - 0.000173 * T) * Csalt)
    return RSwat