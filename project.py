import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# moments of inertia
# symmetric moments
I = 0.02
# moment along symmetry axis
J = 0.01

# parameters for gravitational potential energy
m = 0.5
g = 9.81
ell = 0.1

# gravitational potential energy if theta = 0
U_0 = m*g*ell

# Initialize all Euler angles except psi_dot_0 
theta_0 = np.pi/9
theta_dot_0 = 0

phi_0 = 0
phi_dot_0 = 0

psi_0 =0 

# initialize spin rates
omega_3s = [5, 15, 30]

# set placeholder values for momenta
p_theta_0 = I * theta_dot_0
p_phi_0 = 0
p_psi_0 = 0


def derive_p_theta(theta):
    Num_1 = p_phi_0-p_psi_0*np.cos(theta)
    Num_2 = p_phi_0*np.cos(theta)-p_psi_0
    Denom = I*(np.sin(theta))**3
    T1 = Num_1*Num_2/Denom
    T2 = U_0*np.sin(theta)
    return T1+T2

def derive_phi(theta):
    Num = p_phi_0-p_psi_0*np.cos(theta)
    Denom = I*np.sin(theta)
    T = Num/Denom
    return T

def derive_psi(theta):
    Num_1 = p_psi_0*np.cos(theta)-p_phi_0
    Num_2 = np.cos(theta)
    Denom = I*np.sin(theta)
    T1 = Num_1*Num_2/Denom
    T2 = p_psi_0/J
    return T1+T2


# first order symplectic
def SymplecticEuler(theta_0, ptheta_0, phi_0, psi_0, dt):

    k1 = derive_p_theta(theta_0)
    k2 = derive_phi(theta_0)
    k3 = derive_psi(theta_0)


    ptheta_1 = ptheta_0+k1*dt
    phi_1 = phi_0+k2*dt
    psi_1 = psi_0+k3*dt

    k4 = ptheta_1/I
    theta_1 = theta_0+k4*dt

    return theta_1, ptheta_1, phi_1, psi_1

# first order non-symplectic
def Euler(theta_0, ptheta_0, phi_0, psi_0, dt):

    k1 = derive_p_theta(theta_0)
    k2 = derive_phi(theta_0)
    k3 = derive_psi(theta_0)
    k4 = ptheta_0/I


    ptheta_1 = ptheta_0+k1*dt
    phi_1 = phi_0+k2*dt
    psi_1 = psi_0+k3*dt
    theta_1 = theta_0+k4*dt

    return theta_1, ptheta_1, phi_1, psi_1

# second order symplectic
def Verlet(theta_0, ptheta_0, phi_0, psi_0, dt):
    k1 = derive_p_theta(theta_0)
    k2 = derive_phi(theta_0)
    k3 = derive_psi(theta_0)

    ptheta_1 = ptheta_0 + 0.5*dt*k1
    phi_1 = phi_0 + 0.5*dt*k2
    psi_1 = psi_0 + 0.5*dt*k3

    k4 = ptheta_1/I
    theta_1 = theta_0+k4*dt

    k1 = derive_p_theta(theta_1)
    k2 = derive_phi(theta_1)
    k3 = derive_psi(theta_1)

    ptheta_1 = ptheta_1 + 0.5*dt*k1
    phi_1 = phi_1 + 0.5*dt*k2
    psi_1 = psi_1 + 0.5*dt*k3

    return theta_1, ptheta_1, phi_1, psi_1

# second order non-symplectic
def Midpoint(theta_0, ptheta_0, phi_0, psi_0, dt):

    k11 = derive_p_theta(theta_0)
    k14 = ptheta_0/I

    k21 = derive_p_theta(theta_0+k14*dt/2)
    k22 = derive_phi(theta_0+k14*dt/2)
    k23 = derive_psi(theta_0+k14*dt/2)
    k24 = (ptheta_0+k11*dt/2)/I
    
    ptheta_1 = ptheta_0+dt*k21
    phi_1 = phi_0+dt*k22
    psi_1 = psi_0+dt*k23
    theta_1 = theta_0+dt*k24

    return theta_1, ptheta_1, phi_1, psi_1

# fourth order symplectic
def Yoshida(theta_0, ptheta_0, phi_0, psi_0, dt):
    w0 = 1 / (2 - 2**(1/3))
    w1 = 1 - 2*w0
    theta_1, ptheta_1, phi_1, psi_1 = Verlet(theta_0, ptheta_0, phi_0, psi_0, w0*dt)
    theta_1, ptheta_1, phi_1, psi_1 = Verlet(theta_1, ptheta_1, phi_1, psi_1, w1*dt)
    theta_1, ptheta_1, phi_1, psi_1 = Verlet(theta_1, ptheta_1, phi_1, psi_1, w0*dt)
    return theta_1, ptheta_1, phi_1, psi_1

# fourth order non-symplectic
def RK4(theta_0, ptheta_0, phi_0, psi_0, dt):

    k11 = derive_p_theta(theta_0)
    k12 = derive_phi(theta_0)
    k13 = derive_psi(theta_0)
    k14 = ptheta_0/I

    k21 = derive_p_theta(theta_0+k14*dt/2)
    k22 = derive_phi(theta_0+k14*dt/2)
    k23 = derive_psi(theta_0+k14*dt/2)
    k24 = (ptheta_0+k11*dt/2)/I

    k31 = derive_p_theta(theta_0+k24*dt/2)
    k32 = derive_phi(theta_0+k24*dt/2)
    k33 = derive_psi(theta_0+k24*dt/2)
    k34 = (ptheta_0+k21*dt/2)/I

    k41 = derive_p_theta(theta_0+k34*dt)
    k42 = derive_phi(theta_0+k34*dt)
    k43 = derive_psi(theta_0+k34*dt)
    k44 = (ptheta_0+k31*dt)/I

    
    ptheta_1 = ptheta_0+(dt/6)*(k11+2*k21+2*k31+k41)
    phi_1 = phi_0+(dt/6)*(k12+2*k22+2*k32+k42)
    psi_1 = psi_0+(dt/6)*(k13+2*k23+2*k33+k43)
    theta_1 = theta_0+(dt/6)*(k14+2*k24+2*k34+k44)

    return theta_1, ptheta_1, phi_1, psi_1

for omega_3 in omega_3s:
    psi_dot_0 = omega_3 - phi_dot_0*np.cos(theta_0)

    # set actual parameters for conerved generalized momenta
    p_phi_0 = (I*(np.sin(theta_0))**2+J*(np.cos(theta_0))**2)*phi_dot_0+J*np.cos(theta_0)*psi_dot_0 
    p_psi_0 = J*omega_3

    # time step
    dt = 0.001
    # number of iterations
    N = 1000000
    theta = theta_0
    p_theta = p_theta_0
    theta_dot = theta_dot_0

    phi = phi_0
    p_phi = p_phi_0
    phi_dot = phi_dot_0

    psi = psi_0
    p_psi = p_psi_0
    psi_dot = psi_dot_0

    thetas = [theta]
    p_thetas = [p_theta]
    theta_dots = [theta_dot]

    phis = [phi]
    p_phis = [p_phi]
    phi_dots = [phi_dot]

    psis = [psi]
    p_psis = [p_psi]
    psi_dots = [psi_dot]

    for i in range(N):
        theta, p_theta, phi, psi = SymplecticEuler(theta, p_theta, phi, psi, dt)

        thetas.append(theta)
        p_thetas.append(p_theta)
        theta_dots.append(p_theta/I)


        phi_dot = (1/(I*(np.sin(theta))**2))*(p_phi-p_psi*np.cos(theta))
        phis.append(phi)
        p_phis.append(p_phi)
        phi_dots.append(phi_dot)

        psi_dot = (1/(I*J*(np.sin(theta))**2))*(-p_phi*J*np.cos(theta)+p_psi*(I*(np.sin(theta))**2+J*(np.cos(theta))**2))
        psis.append(psi)
        p_psis.append(p_psi)
        psi_dots.append(psi_dot)

    data = {
    "Theta": thetas,
    "P_Theta": p_thetas,
    "Theta_dot":theta_dots,
    "Phi": phis,
    "P_Phi": p_phis,
    "Phi_dot":phi_dots,    
    "Psi": psis,
    "P_Psi": p_psis,
    "Psi_dot":psi_dots,
    }

    df = pd.DataFrame(data)
    df.to_csv("Project Results/SymplecticEuler/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv", index=False)