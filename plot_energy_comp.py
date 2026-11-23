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


omega_3 =5
N = 10000
dt = 0.001

omega_3 = 30
Ns = [100, 1000, 10000, 100000, 1000000]
Method_Results = ['SymplecticEuler', 'Euler', 'Verlet', 'Midpoint', 'Yoshida', 'RK4']

Emaxss = []

E_0 = (J/2)*omega_3**2+U_0*np.cos(np.pi/9)
for Method_Result in Method_Results:
    Emaxs = []
    for N in Ns:
        df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")

        U = U_0*np.cos(df['Theta'])
        T1 = (I/2)*(df['Theta_dot']**2+df['Phi_dot']**2*(np.sin(df['Theta']))**2)
        T2 = (J/2)*(df['Psi_dot']+df['Phi_dot']*np.cos(df['Theta']))**2
        E = T1+T2+U

        ts = np.linspace(0, int(N*dt), N+1)
        Emaxs.append(np.max(np.abs((E-E_0)/E_0)))
    Emaxss.append(Emaxs)
colors = ['greenyellow', 'olivedrab', 'orchid', 'darkmagenta', 'aqua', 'dodgerblue']
i = 0
for Method_Result in Method_Results:
    plt.plot(np.array(Ns)*dt, Emaxss[i], label=Method_Result+" Method", color = colors[i], marker='o')
    i = i + 1

plt.title('Numerical Method Energy Conservation Comparison')
plt.ylabel('Maximum Normalized Energy Deviation, [(E-E\u2080)/E\u2080]\u2098\u2090\u2093')
plt.xlabel('Simulation Time (s)')
plt.legend()
plt.xscale('log')
plt.yscale('log')
plt.tight_layout()
plt.savefig("Project Results/EnergyChange_Comp.pdf")