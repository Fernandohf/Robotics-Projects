
# Funções e Bibliotecas Utilizadas
from sympy import symbols, pprint, simplify, Eq, diff
from sympy.physics.mechanics import *
from sympy.physics.mechanics.functions import inertia
init_vprinting()

# Variáveis Simbólicas
theta = dynamicsymbols('theta')
dtheta = dynamicsymbols('theta', 1)
omega = dynamicsymbols('omega')
domega = dynamicsymbols('omega', 1)
tau = symbols('tau')
l = symbols('l', positive = True)
r = symbols('r', positive = True)
m, g = symbols('m g')
I_xx, I_yy, I_zz = symbols('I_{xx}, I_{yy}, I_{zz}') 

# Referenciais 
B0 = ReferenceFrame('B0')                         # Referencial Inercial
B1 = B0.orientnew('B1', 'Axis', [theta, B0.z])    # Referencial móvel: theta_1 em relação a B0.z 

# Pontos e Centros de Massa
O = Point('O')                                    # Origem
O.set_vel(B0, 0)              
CM = Point('CM')                                  # Centro de Massa
CM.set_pos(O, r * B1.x)
CM.v2pt_theory(O, B0, B1)

# Corpos Rígidos
I = inertia(B1, I_xx, I_yy, I_zz)                 # Tensor de Inércia
E = RigidBody('E', CM, B1, m, (I, O))             # Elo 1

# Energia Potencial
P = - m * g * B0.y
r_CM = (r * B1.x).express(B0)
E.potential_energy = r_CM.dot(P)

# Forças/Momentos Generalizados
FL = [(B1, tau * B0.z)]
# Método de Lagrange
L = Lagrangian(B0, E).simplify()
LM = LagrangesMethod(L, [theta], frame=B0, forcelist = FL)
L_eq = LM.form_lagranges_equations()
pprint(L_eq)
rhs = LM.rhs()
pprint(rhs)
