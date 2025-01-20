import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from numpy import cos, sin

# System parameters
G = 9.8  # acceleration due to gravity, in m/s^2
L1 = 1.0  # length of pendulum 1 in m
L2 = 0.75  # length of pendulum 2 in m
L3 = 0.5  # length of pendulum 3 in m
L = L1 + L2 + L3  # maximal length of the combined pendulum
M1 = 5.75  # mass of pendulum 1 in kg
M2 = 2.0  # mass of pendulum 2 in kg
M3 = 1.75  # mass of pendulum 3 in kg
t_stop = 15  # how many seconds to simulate
history_len = 500  # how many trajectory points to display

# Rotational damping coefficients for each pendulum
b1 = 0.1
b2 = 0.1
b3 = 0.1

# Initial conditions
th1 = 30.0
w1 = 0.0
th2 = -10.0
w2 = 0.0
th3 = 10.0
w3 = 0.0

state = np.radians([th1, w1, th2, w2, th3, w3])

# Time array
dt = 0.01
t = np.arange(0, t_stop, dt)

# Function for computing derivatives
def derivs(t, state):
    dydx = np.zeros_like(state)

    dydx[0] = state[1]

    delta = state[2] - state[0]
    den1 = (M1 + M2) * L1 - M2 * L1 * cos(delta) * cos(delta)
    dydx[1] = ((M2 * L1 * state[1]**2 * sin(delta) * cos(delta)
                + M2 * G * sin(state[2]) * cos(delta)
                + M2 * L2 * state[3]**2 * sin(delta)
                - (M1 + M2) * G * sin(state[0]))
               / den1) - b1 * state[1]  # Damping term added

    dydx[2] = state[3]

    den2 = (L2 / L1) * den1
    dydx[3] = ((-M2 * L2 * state[3]**2 * sin(delta) * cos(delta)
                + (M1 + M2) * G * sin(state[0]) * cos(delta)
                - (M1 + M2) * L1 * state[1]**2 * sin(delta)
                - (M1 + M2) * G * sin(state[2]))
               / den2) - b2 * state[3]  # Damping term added

    dydx[4] = state[5]

    delta2 = state[4] - state[2]
    den3 = (M2 + M3) * L2 - M3 * L2 * cos(delta2) * cos(delta2)
    dydx[5] = ((M3 * L2 * state[3]**2 * sin(delta2) * cos(delta2)
                + M3 * G * sin(state[4]) * cos(delta2)
                - (M2 + M3) * G * sin(state[2]))
               / den3) - b3 * state[5]  # Damping term added

    return dydx

# Integrating the equations
y = np.empty((len(t), 6))
y[0] = state
for i in range(1, len(t)):
    y[i] = y[i - 1] + derivs(t[i - 1], y[i - 1]) * dt

# Compute x and y positions of the pendulums
x1 = L1 * sin(y[:, 0])
y1 = -L1 * cos(y[:, 0])

x2 = L2 * sin(y[:, 2]) + x1
y2 = -L2 * cos(y[:, 2]) + y1

x3 = L3 * sin(y[:, 4]) + x2
y3 = -L3 * cos(y[:, 4]) + y2

# Plot setup
fig = plt.figure(figsize=(5, 4))
ax = fig.add_subplot(autoscale_on=False, xlim=(-L, L), ylim=(-L, 1.))
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

# Animation function with fading effect for the trace
def animate(i):
    thisx = [0, x1[i], x2[i], x3[i]]
    thisy = [0, y1[i], y2[i], y3[i]]

    history_x = x3[:i]
    history_y = y3[:i]

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)

    # Gradually fade the trace
    alpha = 1 - (i / len(t))
    trace.set_alpha(alpha)

    time_text.set_text(time_template % (i * dt))
    return line, trace, time_text

# Create animation
ani = animation.FuncAnimation(fig, animate, len(y), interval=dt * 1000, blit=True)
#ani.save('/Users/siddharthsharma/myenv312/dynamic_sim_project/tripendulum_animation.mp4', writer='ffmpeg', fps=30)

plt.show()