import sys

print("PYTHON VERSION =", sys.version)
print("sys.path =", sys.path)

import sympy
print("Sympy version =", sympy.__version__)
print("Sympy file =", sympy.__file__)

from sympy import symbols, dynamicsymbols, Function, lambdify
from sympy.physics.mechanics import *
import numpy as np

# Number of pendulums
n = 3

# Generalized coordinates and speeds
q = dynamicsymbols('q:' + str(n + 1))  # e.g. q0, q1, q2, q3
u = dynamicsymbols('u:' + str(n + 1))  # e.g. u0, u1, u2, u3
f = dynamicsymbols('f')                # External force (cart force)

# Define mass and length symbols
m = symbols('m:' + str(n + 1))         # m0, m1, m2, m3
l = symbols('l:' + str(n))             # l0, l1, l2
g, t = symbols('g t')

# Reference frame, points, etc.
I = ReferenceFrame('I')
O = Point('O')
O.set_vel(I, 0)

P0 = Point('P0')
P0.set_pos(O, q[0]*I.x)
P0.set_vel(I, u[0]*I.x)
Pa0 = Particle('Pa0', P0, m[0])

frames = [I]
points = [P0]
particles = [Pa0]
forces = [
    (P0, f*I.x - m[0]*g*I.y)  # force on the cart mass
]
kindiffs = [q[0].diff(t) - u[0]]

# Build out frames/points/particles/forces for each pendulum
for i in range(n):
    Bi = I.orientnew('B' + str(i), 'Axis', [q[i+1], I.z])
    Bi.set_ang_vel(I, u[i+1]*I.z)
    frames.append(Bi)

    Pi = points[-1].locatenew('P' + str(i+1), l[i]*Bi.x)
    Pi.v2pt_theory(points[-1], I, Bi)
    points.append(Pi)

    Pai = Particle('Pa' + str(i+1), Pi, m[i+1])
    particles.append(Pai)

    # Gravity + same external force f * I.x (this is one possible assumption)
    gravity_force = -m[i+1]*g*I.y
    input_force   = f*I.x
    f_i = (Pi, gravity_force + input_force)
    forces.append(f_i)

    # Kinematic ODE:  q'[i+1] - u[i+1] = 0
    kindiffs.append(q[i+1].diff(t) - u[i+1])

# Use Kane's Method to form the equations
kane = KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kindiffs)
fr, frstar = kane.kanes_equations(forces, particles)

# mass_matrix and forcing vector
MM = kane.mass_matrix  # Symbolic mass matrix
FF = kane.forcing      # Symbolic forcing column vector