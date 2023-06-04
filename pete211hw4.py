import math

#Constants
H0 = 80 #m. Piezometric head directly behind pump.
H1 = 55 #m. Piezometric head in the second case.
D = 0.4 #m. Pipe diameter.
R = 0.4 / 2 #m. Pipe radius.
AREA = math.pi * R**2 #m2. Pipe cross sectional area.
EC = 1 * 10**(-3) / D #Relative roughness for concrete pipe.
EP = 0.45 * 10**(-4) / D #Relative roughness for plastic pipe.
FFC = 0.024 #Friction Factor for concrete pipe at EC.
FFP = 0.012 #Friction Factor for plasticpipe at EC.
Q0 = 250 * 10**(-3) #m3/s. Flow rate at first.
Q1 = Q0 * 1.5 #m3/s. Flow rate in the second case.
UW = 1 * 10**(-6) #m2/s. Kinematic viscosity for water.
UO = 2 * 10**(-4) #m2/s. Kinematic viscosity for oil.
L = 2.5 * 10**(3) #m. Pipe lenght.
G = 9.81 #m/s2. Gravitational acceleration

#Velocity
def velocity(qx):
    V = qx / AREA
    print(V)

    return V

#Reynold's Number
def reynold(qx, ux):
    Re = velocity(qx) * D / ux
    print(reynold(Q0, UW))

    return Re

#Question 2A
def question2a():
    hf = FFC * (L/D) * ((velocity(Q0)**2) / (2 * G))
    h = H0 - hf
    print(hf, EC)

    return h

#Question 2B
def question2b():
    hf = FFP * (L/D) * ((velocity(Q0)**2) / (2 * G))
    h = H0 - hf
    print(hf, EP)

    return h

#Question 3
def question3():
    hf = FFC * (L/D) * ((velocity(Q1)**2) / (2 * G))
    h = H0 - hf
    print(Q1)
    print(hf, EC)

    return h

#Question 4
def question4():
    veqn = (H0 - H1) * D * 2 * G / L
    Reqn = D / UW
    
    #Fully Turbulent Flow Assumption
    fassume = 1.325 / (math.log(EP / (3.7 * D)))**2
    vassume = math.sqrt(veqn / fassume)
    reassume = Reqn * vassume
    print(vassume, reassume)

    fassume2 = 0.013 #Needs manual tweaking. Otherwise will give incorrect results.
    vassume2 = math.sqrt(veqn / fassume2)
    reassume2 = Reqn * vassume2
    print(vassume2, reassume2)

    fassume3 = 0.0135 #Needs manual tweaking. Otherwise will give incorrect results.
    vassume3 = math.sqrt(veqn / fassume3)
    reassume3 = Reqn * vassume3
    print(vassume3, reassume3)

    hfcontrol = fassume3 * L / D * (vassume3**2 / (2 * G))
    print(hfcontrol)

    q3 = vassume3 * AREA * 1000
    return q3

#Question 5
def question5():
    veqn = (H0 - H1) * D * 2 * G / L
    Reqn = D / UO
    print(veqn, Reqn)
    
    #Fully Turbulent Flow Assumption
    fassume = 1.325 / (math.log(EP / (3.7 * D)))**2
    vassume = math.sqrt(veqn / fassume)
    reassume = Reqn * vassume
    print(fassume, vassume, reassume)

    fassume2 = 0.038 #Needs manual tweaking. Otherwise will give incorrect results.
    vassume2 = math.sqrt(veqn / fassume2)
    reassume2 = Reqn * vassume2
    print(vassume2, reassume2)

    fassume3 = 0.042 #Needs manual tweaking. Otherwise will give incorrect results.
    vassume3 = math.sqrt(veqn / fassume3)
    reassume3 = Reqn * vassume3
    print(vassume3, reassume3)

    hfcontrol = fassume3 * L / D * (vassume3**2 / (2 * G))
    print(hfcontrol)

    q3 = vassume3 * AREA * 1000
    return q3

#Numbers for question. Made for convenience.
q1 = None
q2a = question2a()
q2b = question2b()
q3 = question3()
q4 = question4()
q5 = question5()

print(q1)