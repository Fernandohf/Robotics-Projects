
# Funções das Bibliotecas Utilizadas
from sympy import symbols, trigsimp, latex, pprint
from sympy.physics.mechanics import dynamicsymbols
from sympy.physics.vector import ReferenceFrame, Vector
from sympy.physics.vector import time_derivative

# Variáveis Simbólicas
theta_1, theta_2, theta_3 = dynamicsymbols('theta_1 theta_2 theta_3')
l_1, l_2 = symbols('l_1 l_2', positive = True)

# Referenciais 
B0 = ReferenceFrame('B0')                 # Referencial Parado
B1 = ReferenceFrame('B1')     
B1.orient(B0, 'Axis', [theta_1, B0.y])    # Referencial móvel: theta_1 em relação a B0.y 
B2 = ReferenceFrame('B2')
B2.orient(B1, 'Axis', [theta_2, B1.z])    # Referencial móvel: theta_2 em relação a B1.z 
B3 = ReferenceFrame('B3')
B3.orient(B2, 'Axis', [theta_3, B2.z])    # Referencial móvel: theta_3 em relação a B2.z 

# Vetores Posição entre os Pontos
B0_r_OA = Vector(0)      # Vetor Nulo
B2_r_AB = l_1 * B2.x     # Vetor que liga os pontos A e B expresso no referencial móvel B2
B3_r_BC = l_2 * B3.x     # Vetor que liga os pontos B e C expresso no referencial móvel B3

# Cinemática do ponto A em relação ao referencial B0
r_A = B0_r_OA
v_A = time_derivative(r_A, B0)
a_A = time_derivative(v_A, B0)

# Cinemática do ponto B em relação ao referencial B0
r_B  = r_A + B2_r_AB.express(B0)
v_B = time_derivative(r_B, B0)
a_B = time_derivative(v_B, B0)

# Cinemática do ponto C em relação ao referencial B0
r_C  = B3_r_BC.express(B0)
v_C = (time_derivative(r_C, B0))
a_C = (time_derivative(v_C, B0))

# Simplifcação dos Resultados
r_A = (r_A.to_matrix(B0)).applyfunc(trigsimp)
v_A = (v_A.to_matrix(B0)).applyfunc(trigsimp)
a_A = (a_A.to_matrix(B0)).applyfunc(trigsimp)
r_B = (r_B.to_matrix(B0)).applyfunc(trigsimp)
v_B = (v_B.to_matrix(B0)).applyfunc(trigsimp)
a_B = (a_B.to_matrix(B0)).applyfunc(trigsimp)
r_C = (r_C.to_matrix(B0)).applyfunc(trigsimp)
v_C = (v_C.to_matrix(B0)).applyfunc(trigsimp)
a_C = (a_C.to_matrix(B0)).applyfunc(trigsimp)

# Resultados de C
pprint(r_C)
pprint(v_C)
pprint(a_C)
