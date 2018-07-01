"""
Adaptador dos códigos em:
http://matplotlib.sourceforge.net/examples/animation/double_pendulum_animated.py
https://jakevdp.github.io/blog/2012/08/18/matplotlib-animation-tutorial/


Author: Fernando Henrique Fernandes
email: nandohfernandes@gmail.com
"""

from numpy import sin, cos, pi
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint, solve_ivp
from scipy.misc import derivative
import matplotlib.animation as animation

class Manipulador1GdL:
    """Classe do Manipulador de 1 Grau de Liberdade

    * init_state: o estado inicial do manipulador [theta, omega] em graus,
        sendo theta e omega são a posição angular e velocidade angular do 
        do braço manipulador.
    * theta_function: o estado alvo [theta, omega] que o controlador deve
        levar o manipulador, sendo theta a  funções que retornam 
        a posição angular em radianos.
    """
    def __init__(self,
                 init_state = [0, 0],
                 theta_function = lambda: 0,  # Função que define a posição angular em função do tempo
                 L = 1.0,                     # Compriemnto do braço (m)
                 R = 0.5,                     # Distância até a Localização do C.M. do braço (m)
                 I = 0.12,                    # Momento de inércia do braço (kg.m^2)
                 M = 1.0,                     # Massa do braço (kg)
                 G = 9.8,                     # Aceleração da gravidade (m/s^2)
                 T = 0.0,                     # Torque aplicado
                 lamb = 25,                   # Valor lambda da equação diferencial do erro 
                 sat = False,                 # Caso exista saturação no motor
                 sat_lim = [-8, 8],           # Limites de Saturação
                 friction = False,             # Caso considere o atrito no modelo
                 edo_method = 0,              # Método de resolução da EDO (LSODA=0, RK45=1)
                                              #     OBS.: LSODA é consideravelmente mais rápido,
                                              #     principalmente no caso com atrito.
                 origin = (0, 0)):
        self.edo_method = edo_method
        self.target_state = [0, 0]
        self.sat = sat
        self.sat_limit = sat_lim
        self.friction = friction
        self.lambda_ = lamb
        self.theta_function = theta_function
        self.init_state = np.asarray(init_state, dtype='float')
        self.params = (L, R, I, M, G)
        self.torque = T
        self.origin = origin
        self.time_elapsed = 0
        self.state = np.deg2rad(self.init_state)    # Converte os ângulos em radianos
        self.error = [0, 0]                         # Lista com os erros [theta_hat, omega_hat]
    
    # Função utilizada pelo solver da EDO
    def f_ode(self, t, y):
        theta, omega = y
        (_, R, I, M, G) = self.params
        T = self.torque
        # Lista com a conversão da equação diferencial de 2 ordem por 2 de primeira ordem
        if self.friction:
            derivs = [omega,
                      1 / I * (T - M * G * R * sin(theta) - 0.25 * np.sign(omega))]
        else:
            derivs = [omega,
                      1 / I * (T - M * G * R * sin(theta))]
        return derivs

    def get_params(self):
        """ Retorna os parâmetros atuais sendo usando no manipulador"""
        return self.params

    def set_params(self, parameters):
        """ Atualiza os parâmetros atuais sendo usando no manipulador"""
        self.params = parameters

    def set_torque(self, t):
            """ Atualiza o parâmetro T (Torque) do manipulador, essa é a 
            maneira que o controlador atua no manipulador, também limita 
            o valor que pode ser definido"""
            if self.sat:
                self.torque = max(-5.0, min(t, 5.0))
            else:
                self.torque = t

    def position(self):
        """ Retorna a posição (x,y) atual da ponta do manipulador"""
        L = self.params[0]
        x = np.cumsum([self.origin[0], L * cos(self.state[0])])
        y = np.cumsum([self.origin[1], L * sin(self.state[0])])
        return (x, y)
    
    def target_position(self):
            """ Retorna a posição alvo (x,y) atual da ponta do manipulador"""
            L = self.params[0]
            x = np.cumsum([self.origin[0], L * cos(self.target_state[0])])
            y = np.cumsum([self.origin[1], L * sin(self.target_state[0])])
            return (x, y)

    def control_arm(self):
        """Função que calcular o torque necessário para alcaçar os estados desejados """
        # Valores atuais
        theta_curr, omega_curr = self.state

        # Tranforma as funções alvo em valores alvo
        theta_target = self.theta_function(self.time_elapsed)
        omega_target = derivative(self.theta_function, self.time_elapsed, 1e-18)
        alfa_target = derivative(self.theta_function, self.time_elapsed, 1e-18, 2)
        self.target_state[0] = theta_target
        self.target_state[1] = omega_target

        # Calcula os erros
        theta_hat = theta_curr - theta_target
        omega_hat = omega_curr - omega_target
        self.error[0] = theta_hat
        self.error[1] = omega_hat
        
        # Define o Torque necessário
        (_, R, I, M, G) = self.params
        lambda_ = self.lambda_
        torq = (M * G * R * sin(theta_curr) +
               I * (alfa_target - 2 * lambda_ * omega_hat - (lambda_ ** 2) * theta_hat))
        self.set_torque(torq)

    def step(self, dt):
        """Executa um passo de comprimento dt, invoca control _arm e 
        atualiza o estado atual do manipulador"""
        # Atualiza o torque das EDO's
        self.control_arm()
        
        # Tempos 
        t0 = self.time_elapsed
        tf = t0 + dt
        
        if self.edo_method == 0:
            # Solucionador padrão de EDO do Scipy - LSODA
            self.state = odeint(self.f_ode, self.state, [t0, tf], tfirst=True)[1]
        else:
            # Runge-Kutta 4 ordem
            kr45 = solve_ivp(self.f_ode, (t0, tf),
                             self.state, t_eval=[tf])
            self.state[0] = kr45.y[0][0]
            self.state[1] = kr45.y[1][0]
        #print(self.state)
        self.time_elapsed += dt


#------FUNÇÃO POSIÇÃO ANGULAR E INICIALIZAÇÃO DA CLASSE------#
# Função que retorna a posição angular (em radianos) em função do tempo
def f_t(t):
    #angulo = np.deg2rad(30)
    angulo = pi / 6 * (1 - cos(2 * pi * t))
    # if t < 2:
    #     angulo = np.deg2rad(30)
    # elif t < 6:
    #     angulo = np.deg2rad(180)
    # elif t < 10:
    #     angulo = np.deg2rad(-20)
    # else:
    #     angulo = np.deg2rad(-90)
    return angulo

# Inicializa a classe
manipulador = Manipulador1GdL([0, 0], f_t, edo_method=1, friction=True)
dt = 1./30              # Intevalor de tempo  = 30 fps
animate_model = True    # Se plotará a animação, caso contrário, apenas os gráficos
tmin, tmax = [0 , 10.0] # Tempo mínimo e máximo
torques = []            # Torques obtidos a cada iteração
thetas = []             # Estados dos thetas a cada iteração
thetas_hat = []         # Erros do thetas obtidos a cada iteração
omegas = []             # Estados dos omegas a cada iteração
omegas_hat = []         # Erros do omega obtidos a cada iteração
times = []              # Tempo de cada iteração

#------ GRÁFICOS E ANIMAÇÃO------#
# Figura
fig = plt.figure(figsize=(8, 8))

# Gráficos
if animate_model:
    # Posição dos Subplots com animação
    grid = plt.GridSpec(2, 4, wspace=0.2, hspace=0.6)
    pos_anim = grid[0:2, 0:2]
    pos_torq = grid[0, 2:]
    pos_error = grid[1, 2:]
    sim_ax = fig.add_subplot(pos_anim, aspect='equal', autoscale_on=False,
                             xlim=(-1.5, 1.5), ylim=(-1.5, 1.5))
    sim_ax.grid(alpha=0.8)
    sim_ax.set_title("Animação do Manipulador")

    torq_ax = fig.add_subplot(pos_torq, autoscale_on=False,
                            xlim=(0, 6), ylim=(manipulador.sat_limit[0], manipulador.sat_limit[1]))
    torq_ax.set_title("Torque x Tempo")
    torq_ax.get_xaxis().set_visible(False)

    error_ax = fig.add_subplot(pos_error, autoscale_on=False,
                            xlim=(0, 6), ylim=(-5, 5))
    error_ax.set_title(r"[$\hat{\theta},\dot{\hat{\theta}}$] x Tempo")
    error_ax.get_xaxis().set_visible(False)
    
    # Objetos Desenhados
    line_anim, = sim_ax.plot([], [], 'bo-', lw=2)
    line_ideal, = sim_ax.plot([], [], 'co-', lw=2)
    line_torq, = torq_ax.plot([], [], 'r-', lw=2)
    line_error_theta, = error_ax.plot([], [], 'g-', lw=2, label=r"$\hat{\theta}$")
    line_error_omega, = error_ax.plot([], [], 'y-', lw=2, label=r"$\dot{\hat{\theta}}$")
    time_text = sim_ax.text(0.04, 0.04, '', transform=sim_ax.transAxes)
    error_ax.legend()

else:
    grid = plt.GridSpec(4, 2, wspace=0.4, hspace=0.5)
    pos_torq = grid[0, 0:]
    pos_error = grid[1, 0:]
    pos_error2 = grid[2, 0:]
    pos_error3 = grid[3, 0:]
    
    torq_ax = fig.add_subplot(pos_torq)
    torq_ax.set_title("Torque x Tempo")

    error_ax = fig.add_subplot(pos_error)
    error_ax.set_title(r"[$\theta,\hat{\theta}$] x Tempo")

    error2_ax = fig.add_subplot(pos_error2)
    error2_ax.set_title(r"[$\dot{\theta},\dot{\hat{\theta}}$] x Tempo")

    error3_ax = fig.add_subplot(pos_error3)
    error3_ax.set_title(r"$\dot{\hat{\theta}}$ x $\hat{\theta}$")


def init():
    """initialize animation"""
    line_anim.set_data([], [])
    line_ideal.set_data([], [])
    line_torq.set_data([], [])
    line_error_theta.set_data([], [])
    line_error_omega.set_data([], [])
    time_text.set_text('')

    return line_anim, line_ideal, line_torq, line_error_theta, line_error_omega, time_text


def animate(i):
    """perform animation step"""
    global manipulador, dt, torques, thetas_hat, omegas_hat, times, torq_ax, error_ax

    manipulador.step(dt)

    # Atualiza os dados
    t_ = manipulador.time_elapsed
    thetas_hat.append(manipulador.error[0])
    omegas_hat.append(manipulador.error[1])
    torques.append(manipulador.torque)
    times.append(t_)
    xmin, xmax = torq_ax.get_xlim()
    
    # Ajusta o eixo do tempo, x
    if t_ > xmax:
        torq_ax.set_xlim((xmin, xmax * 2))
        torq_ax.figure.canvas.draw()
        error_ax.set_xlim((xmin, xmax * 2))
        error_ax.figure.canvas.draw()
    
    # Atualiza os objetos desenhados com os dados obtidos
    line_anim.set_data(*manipulador.position())
    line_ideal.set_data(*manipulador.target_position())
    line_torq.set_data(times, torques)
    line_error_theta.set_data(times, omegas_hat)
    line_error_omega.set_data(times, thetas_hat)
    time_text.set_text('Tempo = %.1f' % manipulador.time_elapsed)
    
    return line_anim, line_ideal, line_torq, line_error_theta, line_error_omega, time_text

if animate_model:
    # Animação e Simulação dos Resultados
    from time import time
    t0 = time()
    animate(0)
    t1 = time()
    interval = 1000 * dt - (t1 - t0)  # Intervalo baseado em dt e no tempo para se calcular um passo

    ani = animation.FuncAnimation(fig, animate, frames=300,
                                interval=interval, blit=True, init_func=init)
else:
    # Gráficos dos Resultados
    time_list = np.arange(tmin, tmax, dt)
    for deltat in time_list:
        manipulador.step(dt)
        thetas.append(manipulador.state[0])
        thetas_hat.append(manipulador.error[0])
        omegas.append(manipulador.state[1])
        omegas_hat.append(manipulador.error[1])
        torques.append(manipulador.torque)
        
    # Desenha os Objetos
    line_torq, = torq_ax.plot(time_list, torques, 'r-', lw=2)
    line_theta, = error_ax.plot(time_list, thetas, 'g-', lw=2, label=r"$\theta$")
    line_error_theta, = error_ax.plot(time_list, thetas_hat, 'y-', lw=2, label=r"$\hat{\theta}$")
    line_omega, = error2_ax.plot(time_list, omegas, 'b-', lw=2, label=r"$\dot{\theta}$")
    line_error_omega, = error2_ax.plot(time_list, omegas_hat, 'c-', lw=2, label=r"$\dot{\hat{\theta}}$")
    line_erros_omega_theta, = error3_ax.plot(omegas_hat, thetas_hat, 'm-', lw=2, label=r"$\dot{\hat{\theta}}$")
    error_ax.legend()
    error2_ax.legend()

plt.show()

