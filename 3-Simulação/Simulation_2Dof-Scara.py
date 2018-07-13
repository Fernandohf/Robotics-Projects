"""
Adaptador dos códigos em:
http://matplotlib.sourceforge.net/examples/animation/double_pendulum_animated.py
https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/


Author: Fernando Henrique Fernandes
email: nandohfernandes@gmail.com
"""

from numpy import sin, cos, pi
import dill, os, time
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from scipy.misc import derivative
import matplotlib.animation as animation
from sympy import symbols, pprint, simplify, Derivative, lambdify, Dummy
from sympy.physics.mechanics import *
from sympy.physics.mechanics.functions import inertia
dill.settings['recurse'] = True

class Manipulador2GdL:
    """Classe do Manipulador de 2 Grau de Liberdade

    * init_state: o estado inicial do manipulador [theta_1, theta_2, omega_1 e omega_2] em graus ou rad/s,
        sendo thetas e omegas as posições angulares e velocidades angulares dos elos do manipulador.
    """
    def __init__(self,
                 init_state = [0, 0, 0, 0],
                 L1 = 1.0,                    # Compriemnto dos braços 1 (m)
                 L2 = 0.6,                    # Compriemnto dos braços 2 (m)
                 R1 = 0.5,                    # Distância até a Localização do C.M. do braço 1(m)
                 R2 = 0.3,                    # Distância até a Localização do C.M. do braço 2(m)
                 I1 = 0.15,                   # Momento de inércia do braço 1(kg.m^2)
                 I2 = 0.05,                   # Momento de inércia do braço 2(kg.m^2)
                 M1 = 1.0,                    # Massa do braço 1(kg)
                 M2 = 0.5,                    # Massa do braço 2(kg)
                 G = 9.81,                    # Aceleração da gravidade (m/s^2)
                 T1 = 0.0,                    # Torque aplicado no elo 1
                 T2 = 0.0,                    # Torque aplicado no elo 2
                 edo_method = 'LSODA',        # Método de resolução da EDO (LSODA=0, RK45=1)
                                              #     OBS.: LSODA é consideravelmente mais rápido,
                                              #     principalmente no caso com atrito.
                 origin = (0, 0)):
        self.edo_method = edo_method
        self.init_state = np.asarray(init_state, dtype='float')
        self.num_params = (L1, L2, R1, R2, I1, I2, M1, M2, G)
        self.torque = [T1, T2]
        self.origin = origin
        # Converte os angulos em radianos e no referencia vertical
        self.thetas_1 = []
        self.thetas_2 = []
        self.omegas_1 = []
        self.omegas_2 = []
        self.time = []        
        self.rhs_lambdified = []
        # frames por segundo da animação
        self.fps = 1./60    
        self.solution_edo = False
        # Calcula a Dinâmica do Corpo
        self.solve_dynamics()

    # Solução da Dinâmica do Manipulador
    def solve_dynamics(self):    
        """Calcula a solução da dinâmica do manipulador e já coloca em
         formato pronto para ser calculado """
        # Variáveis Simbólicas do problema
        theta_1, theta_2 = dynamicsymbols('theta_1 theta_2')
        dtheta_1, dtheta_2 = dynamicsymbols('theta_1 theta_2', 1)
        tau_1, tau_2 = symbols('tau_1 tau_2')
        symb_dynamics = (theta_1, theta_2, dtheta_1, dtheta_2, tau_1, tau_2)

        l_1, l_2 = symbols('l_1 l_2', positive = True)
        r_1, r_2 = symbols('r_1 r_2', positive = True)
        m_1, m_2, g = symbols('m_1 m_2 g')
        I_zz_1, I_zz_2 = symbols('I_{1zz} I_{2zz}') 
        symb_params = (l_1, l_2, r_1, r_2, I_zz_1, I_zz_2, m_1, m_2, g)
        # Nome do arquivo salvo com a equação do movimento
        here = os.path.dirname(os.path.realpath(__file__))
        file_name = (here + "\\models\\scara_dynamics_" + 
                    "".join(str(e) + "_" for e in self.num_params) + ".data")
        try:
            rhs = dill.load(open(file_name, "rb"))
            self.rhs_lambdified = rhs
            return
        except:
            print("Modelo dinâmico ainda não foi calculado para esses parâmetros. Aguarde...")
        initial_time = time.time()

        # Referenciais 
        B0 = ReferenceFrame('B0')                 # Referencial Parado
        B1 = ReferenceFrame('B1')     
        B1.orient(B0, 'Axis', [theta_1, B0.z])    # Referencial móvel: theta_1 em relação a B0.z 
        B2 = ReferenceFrame('B2')
        B2.orient(B1, 'Axis', [theta_2, B1.z])    # Referencial móvel: theta_2 em relação a B1.z 

        # Pontos e Centros de Massa
        O = Point('O')
        O.set_vel(B0, 0)
        A = Point('A')     
        A.set_pos(O, l_1 * B1.x)                     
        A.v2pt_theory(O, B0, B1)
        CM_1 = Point('CM_1')
        CM_1.set_pos(O, r_1 * B1.x)
        CM_1.v2pt_theory(O, B0, B1)
        CM_2 = Point('CM_2')
        CM_2.set_pos(A, r_2 * B2.x)
        CM_2.v2pt_theory(O, B0, B2)

        # Corpos Rígidos
        I_1 = inertia(B1, 0, 0, I_zz_1)            
        E_1 = RigidBody('Elo_1', CM_1, B1, m_1, (I_1, CM_1)) # Elo 1
        I_2 = inertia(B2, 0, 0, I_zz_2)                 
        E_2 = RigidBody('Elo_2', CM_2, B2, m_2, (I_2, CM_2)) # Elo 2) 

        # Energia Potencial
        P_1 = -m_1 * g * B0.y
        r_1_CM = (r_1 * B1.x).express(B0)
        E_1.potential_energy = r_1_CM.dot(P_1)
        P_2 = -m_2 * g * B0.y
        r_2_CM = (r_2 * B2.x).express(B0)
        E_2.potential_energy = r_2_CM.dot(P_2)

        # Forças/Momentos Generalizados
        FL = [(B1, tau_1 * B0.z),(B2, tau_2 * B0.z)]
        # Método de Lagrange
        L = Lagrangian(B0, E_1, E_2).simplify()
        LM = LagrangesMethod(L, [theta_1, theta_2], frame=B0, forcelist=FL)
        motion_eq = LM.form_lagranges_equations()
        rhs = LM.rhs()
        # Salva as Equação em formato fácil de se obter Solução já com os Parâmetros
        dummys = [Dummy() for i in symb_dynamics]
        rhs = msubs(rhs, dict(zip(symb_params, self.num_params)))
        dummydict = dict(zip(symb_dynamics, dummys))
        rhs = msubs(rhs, dummydict)
        
        # Lambdify as equações
        for f in rhs:
            self.rhs_lambdified.append(lambdify(dummys, f))

        # Salva as equações lambidificadas e já com as constantes subtituídas no disco
        os.makedirs(os.path.dirname(file_name), exist_ok=True)
        dill.dump(self.rhs_lambdified, open(file_name, "wb"))

        print("Modelo calculado e salvo em " + file_name + ".Tempo: " +
               str(time.time() - initial_time))

    # Função utilizada pelo solver da EDO
    def f_ode(self, t, y):
        theta_1, theta_2, omega_1, omega_2 = y
        T1, T2 = self.torque
        # Lista com a conversão da equação diferencial do Sympy em função solucionável do Scipy
        return [f(theta_1, theta_2, omega_1, omega_2, T1, T2) for f in self.rhs_lambdified]

    def set_params(self, parameters):
        """ Atualiza os parâmetros atuais sendo usando no manipulador"""
        self.params = parameters

    def position(self, theta1, theta2):
        """ Retorna a posição (x,y) atual da ponta do manipulador"""
        L1 = self.num_params[0]
        L2 = self.num_params[1]
        x = np.cumsum([self.origin[0], L1 * cos(theta1), L2 * cos(theta2)])
        y = np.cumsum([self.origin[1], -L1 * sin(theta1), -L2 * sin(theta2)])
        return (x,y)

    def solve_edo(self, dt = 10):
        """Executa asimulaçao do modelo em um intervalor de tempo dt"""
        # Postos que a função será avaliada
        self.time = np.arange(0, dt, self.fps)
        try:
            # Solucionador padrão de EDO do Scipy
            solution = solve_ivp(self.f_ode, (0, dt), self.init_state, method=self.edo_method, t_eval=self.time)
            self.thetas_1 = solution.y[0]
            self.thetas_2 = solution.y[1]
            self.omegas_1 = solution.y[2]
            self.omegas_2 = solution.y[3]
            self.time = solution.t
            self.solution_edo = solution.success
        except:
            print("A EDO não pôde ser selecionada, verifique se a solução da dinâmica" +
                  "do corpo já foi calculada através do método solve_dynamics().")

    def simulate_model(self):
        """Simula o manipulador a partir das equações dinâmicas"""
        # Verifica se existe solução calculada
        if not self.solution_edo:
            print("As soluções da EDOs devem ser calculadas antes!")
            return

        #------ GRÁFICOS E ANIMAÇÃO ------#
        # Variáveis iniciais
        dt = self.fps       # fps
        thetas_1 = []       # Estados dos thetas a cada iteração
        thetas_2 = []       # Estados dos thetas a cada iteração
        omegas_1 = []       # Estados dos omegas a cada iteração
        omegas_2 = []       # Estados dos omegas a cada iteração
        times = []          # Tempo de cada iteração
        
        # Figura
        fig = plt.figure(figsize=(10, 5))

        # Posição dos Subplots com animação
        time_max = 8
        grid = plt.GridSpec(2, 4, wspace=0.4, hspace=0.4)
        pos_anim = grid[0:, 0:-2]
        pos_state_1 = grid[0, 2:]
        pos_state_2 = grid[1, 2:]

        #Definição dos Subplots em sí
        sim_ax = fig.add_subplot(pos_anim, aspect='equal', autoscale_on=False,
                                xlim=(-2.5, 2.5), ylim=(-2.5, 2.5))
        sim_ax.grid(alpha=0.6)
        sim_ax.set_title("Animação do Manipulador")

        state_ax_1 = fig.add_subplot(pos_state_1, autoscale_on=False,
                                xlim=(0, time_max), ylim=(-5, 5))
        state_ax_1.set_title(r"[$\theta_1, \theta_2$] x Tempo")
        state_ax_1.get_xaxis().set_visible(False)

        state_ax_2 = fig.add_subplot(pos_state_2, autoscale_on=False,
                                xlim=(0, time_max), ylim=(-15, 15))
        state_ax_2.set_title(r"[$\dot{\theta}, \dot{\theta_2}$] x Tempo")
        state_ax_2.get_xaxis().set_visible(False)
            
        # Objetos Desenhados
        line_anim, = sim_ax.plot([], [], 'bo-', lw=2, label="Posição do Manipulador")
        time_text = sim_ax.text(0.04, 0.04, '', transform=sim_ax.transAxes)
        sim_ax.legend()

        line_state_theta_1, = state_ax_1.plot([], [], 'tab:purple', lw=2, label=r"$\theta_1$")
        line_state_theta_2, = state_ax_1.plot([], [], 'tab:pink', lw=2, label=r"$\theta_2$")
        state_ax_1.legend()

        line_state_omega_1, = state_ax_2.plot([], [], 'tab:orange', lw=2, label=r"$\dot{\theta_1}$")
        line_state_omega_2, = state_ax_2.plot([], [], 'tab:red', lw=2, label=r"$\dot{\theta_2}$")
        state_ax_2.legend()

        # Funções necessárias para a animação pelo matplolib
        def init():
            """initialize animation"""
            line_anim.set_data([], [])
            line_state_theta_1.set_data([], [])
            line_state_theta_2.set_data([], [])
            line_state_omega_1.set_data([], [])
            line_state_omega_2.set_data([], [])
            time_text.set_text('')

            return line_anim, line_state_theta_1, line_state_theta_2, line_state_omega_1, line_state_omega_2, time_text
        
        def animate(i):
            """perform animation step"""
            theta_1 = self.thetas_1[i]
            theta_2 = self.thetas_2[i]
            omega_1 = self.omegas_1[i]
            omega_2 = self.omegas_2[i]
            t_ = self.time[i]
            # Atualiza os dados
            thetas_1.append(theta_1)
            thetas_2.append(theta_2)
            omegas_1.append(omega_1)
            omegas_2.append(omega_2)
            times.append(t_)
            xmin, xmax = state_ax_1.get_xlim()
            
            # Ajusta o eixo do tempo, x
            if t_ > xmax:
                state_ax_1.set_xlim((xmin, xmax * 2))
                state_ax_2.set_xlim((xmin, xmax * 2))
                state_ax_1.figure.canvas.draw()
            
            # Atualiza os objetos desenhados com os dados obtidos
            line_anim.set_data(self.position(theta_1, theta_2))
            line_state_theta_1.set_data(times, thetas_1)
            line_state_theta_2.set_data(times, thetas_2)
            line_state_omega_1.set_data(times, omegas_1)
            line_state_omega_2.set_data(times, omegas_2)
            time_text.set_text('Tempo = %.1f' % t_)
            
            return line_anim, line_state_theta_1, line_state_theta_2, line_state_omega_1, line_state_omega_2, time_text

        t0 = time.time()
        animate(0)
        tf = time.time()
        deltat = tf - t0
        print(deltat)
        # Executa a animação
        interval = 1000 * dt - (deltat * 5)     # Intervalo em ms
        fr = len(self.time)
        ani = animation.FuncAnimation(fig, animate, frames=fr, interval=interval,
                                      blit=True, init_func=init)
        plt.show()


#-------- EXECUÇÃO ---------#
# Inicializa a classe
manipulador = Manipulador2GdL([0, 0, 0, 0], edo_method='RK45')
manipulador.solve_edo(20)
manipulador.simulate_model()
