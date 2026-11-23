import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.signal import find_peaks

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

# Numerical parameters
N = 10000
dt = 0.001
Method_Results = ['SymplecticEuler', 'Euler', 'Verlet', 'Midpoint', 'Yoshida', 'RK4']
# select which numerical method result should be used
Method_Result = Method_Results[0]

omega_3s = [5,15,30]

# for part 2


# values that need to be plotted
values = [['Theta', '\u03b8(t)', "Angle (\u00b0)", "Angle (rad)"], 
          ['Phi', '\u03d5(t)', "Angle (\u00b0)", "Angle (rad)"], 
          ['Phi_mod', '\u03d5(t)', "Angle (\u00b0)", "Angle (rad)"],
          ['Phi_dot', '\u03d5\u0307(t)', "Angular Velocity (\u00b0/s)", "Angular Velocity (rad/s)"]]


for value in values:
    fig, ax = plt.subplots(layout='constrained')
    for omega_3 in omega_3s:
        df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")

        ts = np.linspace(0, int(N*dt), N+1)
        angles = []
        if value[0] == 'Phi_mod':
            angles = np.mod(df['Phi']*(180/np.pi),360)
        else:    
            angles = df[value[0]]*(180/np.pi)
        ax.plot(ts, angles, label="\u03C9\u2083="+str(omega_3)+" rad/s")


    ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1))
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.15),
        ncols=3
    )
    plt.xlabel("Time (s)")
    plt.ylabel(value[2])
    secax = ax.secondary_yaxis('right', functions=(np.deg2rad, np.rad2deg))
    secax.set_ylabel(value[3])
    plt.title(value[1]+" Comparison for Differing Spin Rates")

    #save results as a figure 
    plt.savefig("Project Results/"+Method_Result+value[0]+".pdf")

    plt.clf()

# plot theta dot against theta
fig, ax = plt.subplots(layout='constrained')
for omega_3 in omega_3s:
    df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")
    angles = df['Theta']*(180/np.pi)
    angledots = df['Theta_dot']*(180/np.pi)
    ax.plot(angles, angledots, label="\u03C9\u2083="+str(omega_3)+" rad/s")

ax.legend()

plt.xlabel("Angle (\u00b0)")
plt.ylabel("Angular Velocity (\u00b0/s)")
secax = ax.secondary_yaxis('right', functions=(np.deg2rad, np.rad2deg))
secax.set_ylabel('Angular Velocity (rad/s)')
secax1 = ax.secondary_xaxis('top', functions=(np.deg2rad, np.rad2deg))
secax1.set_xlabel('Angle (rad)')
plt.title("Nutation Trajectories (\u03b8, \u03b8\u0307)")

#save results as a figure 
plt.savefig("Project Results/"+Method_Result+"ThetaParametric.pdf")
plt.clf()

# for part 3 
mean_periods = []
theta_maxs = []
theta_mins = []
for omega_3 in omega_3s:
    df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")
    peaks, _ = find_peaks(df['Theta'])
    valleys, _ = find_peaks(-1*df['Theta'])
    theta_max = np.max(df['Theta'][peaks]*(180/np.pi))
    theta_min = np.min(df['Theta'][valleys]*(180/np.pi))
    ts = np.linspace(0, int(N*dt), N+1)
    mean_period = np.average(np.diff(ts[peaks]))    

    mean_periods.append(mean_period)
    theta_maxs.append(theta_max)
    theta_mins.append(theta_min)


#save results as a table
data = {
    "Spin Rate (rad/s)": omega_3s,
    "Mean Period (s)": mean_periods,
    "Max Theta": theta_maxs,
    "Min Theta": theta_mins
    }

df = pd.DataFrame(data)
table_name= "Project Results/Results"+Method_Result+"N"+str(N)+"dt"+str(dt)+".csv"
df.to_csv(table_name, index=False)

# for part 4
fig, ax = plt.subplots(layout='constrained')
maximum_deviations = []
for omega_3 in omega_3s:
    df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")
    E_0 = (J/2)*omega_3**2+U_0*np.cos(np.pi/9)


    U = U_0*np.cos(df['Theta'])
    T1 = (I/2)*(df['Theta_dot']**2+df['Phi_dot']**2*(np.sin(df['Theta']))**2)
    T2 = (J/2)*(df['Psi_dot']+df['Phi_dot']*np.cos(df['Theta']))**2
    E = T1+T2+U

    ts = np.linspace(0, int(N*dt), N+1)
    maximum_deviations.append(np.max(np.abs((E-E_0)/E_0)))
    plt.plot(ts, (E-E_0)/E_0, label="\u03C9\u2083="+str(omega_3)+" rad/s")
plt.title('Change in Energy over Time')
plt.ylabel('Normalized Energy Deviation, (E-E\u2080)/E\u2080')
plt.xlabel('Time (s)')
plt.legend(
    loc="upper center",
    bbox_to_anchor=(0.5, -0.15),
    ncols=3
)

plt.tight_layout()
#save results as a figure 
plt.savefig("Project Results/"+Method_Result+"EnergyChange.pdf")
plt.clf()
#add results to the table
df = pd.read_csv(table_name)
df["Maximum Energy Deviation"] = np.array(maximum_deviations)
df.to_csv(table_name, index=False)

# for part 5
i = 0
colors =['tab:blue', 'tab:orange','tab:green']
for omega_3 in omega_3s:
    df = pd.read_csv("Project Results/"+Method_Result+"/omega"+str(omega_3)+"N"+str(N)+"dt"+str(dt)+".csv")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')


    x = np.sin(df['Theta'])*np.cos(df['Phi'])
    y = np.sin(df['Theta'])*np.sin(df['Phi'])
    z = np.cos(df['Theta'])


    # Plot 3D line and scatter points
    ax.plot(x, y, z, label='Spiral Path', c= colors[i])

    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1.1, 1.1)
    ax.set_zlim(-1.1, 1.1)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    plt.title("Symmetry-Axis Unit Vector (\u03C9\u2083="+str(omega_3)+" rad/s)")

    plt.savefig("Project Results/"+Method_Result+"omega3_"+str(omega_3)+"_3D.pdf")
    plt.clf()

    i = i+1