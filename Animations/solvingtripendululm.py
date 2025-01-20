from sympy import symbols
from sympy.physics.mechanics import (
    dynamicsymbols, ReferenceFrame, Point, Particle, KanesMethod
)

# 1) Number of pendulums
n = 3

# 2) Generalized coordinates & speeds
#    q0 is the cart's horizontal position, q1..q3 are pendulum angles
q = dynamicsymbols('q0:{}'.format(n+1))  # q0, q1, q2, q3
u = dynamicsymbols('u0:{}'.format(n+1))  # u0, u1, u2, u3
f = dynamicsymbols('f')                 # External force on the cart

# 3) Symbolic parameters
m = symbols('m0:{}'.format(n+1))  # m0, m1, m2, m3
l = symbols('l0:{}'.format(n))    # l0, l1, l2
g, t = symbols('g t')

# 4) Inertial frame & origin
I = ReferenceFrame('I')
O = Point('O')
O.set_vel(I, 0)

# 5) Cart point & particle
P0 = Point('P0')
P0.set_pos(O, q[0]*I.x)
P0.set_vel(I, u[0]*I.x)
Pa0 = Particle('Pa0', P0, m[0])  # We still define the cart as a Particle

# 6) Forces list (apply force to the POINT P0, not the Particle)
forces = [
    (P0, f*I.x - m[0]*g*I.y),   # horizontal force f & gravity
]

# 7) Kinematic constraints
kindiffs = [q[0].diff(t) - u[0]]

# Lists to hold points & particles
points = [P0]
particles = [Pa0]

# 8) Build each pendulum link
for i in range(n):
    # Frame for pendulum i
    Bi = I.orientnew('B{}'.format(i), 'Axis', [q[i+1], I.z])
    Bi.set_ang_vel(I, u[i+1]*I.z)
    
    # New point Pi at the end of this link
    Pi = points[-1].locatenew('P{}'.format(i+1), l[i]*Bi.x)
    Pi.v2pt_theory(points[-1], I, Bi)
    points.append(Pi)
    
    # New particle for pendulum bob
    Pai = Particle('Pa{}'.format(i+1), Pi, m[i+1])
    particles.append(Pai)
    
    # Force: gravity + same horizontal input if desired
    gravity_force = -m[i+1]*g*I.y
    input_force   = f*I.x
    # Apply to Pi (the point), not Pai
    forces.append((Pi, gravity_force + input_force))
    
    # Kinematic eq: q'[i+1] - u[i+1] = 0
    kindiffs.append(q[i+1].diff(t) - u[i+1])

# 9) Form equations of motion via Kane’s Method
kane = KanesMethod(I, q_ind=q, u_ind=u, kd_eqs=kindiffs)
fr, frstar = kane.kanes_equations(forces, particles)

MM = kane.mass_matrix
FF = kane.forcing

print("Mass matrix (MM) =\n", MM)
print("\nForcing (FF) =\n", FF)