Fernando Henrique Fernandes

# Problema

Realize o controle da posição do manipulador robótico abaixo usando o método de `Linearização por Realimentação`, sabendo que a sua equação de movimento é dada por $I_{zz}\ddot{\theta} + mgr \text{sen}\theta = \tau$. Para isso, siga os seguintes passos:

![Imagem do Manipulador](https://i.imgur.com/Njacvei.png)

Para isso, siga os seguintes passos:

1. Implementar o método de Runge-Kuttta de 4ª ordem para a simular a dinâmica do manipulador.

2. Projetar um controlador pelo método de linearização por realimentação de modo que a garra vá para a posição $\left(\sqrt{3}/2, 1/2\right)$. 

3. Aproveitando esse mesmo controlador faça com que o manipulador oscile de acordo com a função ${\theta}_d = \frac{\pi}{6}\left(1-\text{cos} 2 \pi t \right)$.

4. Ainda aproveitando o mesmo controlador, verifique a influência que o atrito seco na junta de rotação causa no controle do braço. O atrito é modelo pela adição do termo $0.25 \text{sgn} \theta$ no lado esquerdo da equação de movimento apresentada anteriormente.





# Solução

Toda a implementação necessária para a simulação e resolução do problema proposto foi feita em Python. Para que o modelo ele represente a imagem mostrada no problema, uma pequena alteração foi realizada, pois a equação diferencial foi encontrada considerando que o $\theta$  é medido em relação a parte negativa do eixo $y$, enquanto que na imagem esse valor é medido em relação ao eixo $x$. Então, foi adicionado $\pi/2$ para todos os ângulos fornecidos, de maneira que os ângulo possam ser medidos em relação ao eixo $x$. Além disso, da mesma maneira que o problema foi apresentado, as respostas serão elaboradas por etapas.

## 1. Solução da EDO

Antes da implementação da solução numérica da EDO, temos que ter todas as variáveis que representam as propriedades físicas do modelo. Esses valores foram oferecidos junto com o problema e estão exibidos na tabela abaixo no sistema internacional (SI).

| Variável | Valor      |
|:--------:|:---------- |
| $I_{zz}$ | 0.12 kg m² |
| $m$      | 1 kg       |
| $l$      | 1 m        |
| $r$      | 0.5 m      |
| $g$      | 9.81 m/s²  |

Além disso, também foram fornecidos os valores iniciais  da equação diferencial que são $\theta\left(0\right)=0$ e $\dot{\theta}\left(0\right)=0$. Dessa maneira, podemos iniciar a implementação da solução do método Runge-Kutta de 4ª Ordem. Para isso, procurou-se por alguma biblioteca que realizasse a implementação desse método em Python e algumas foram encontradas. A biblioteca utilizada foi a **Scipy**, que possui vários métodos de soluções numéricas de EDO, inclusive o do Runge-Kutta.

Primeiramente, a equação diferencial ordinária (EDO) de 2ª ordem foi manipulada para que fosse utilizada pelo método. Para tal, é necessário transformá-la em um sistema de EDO de 1ª ordem. Isso é facilmente realizado por manipulação algébrica simples, como mostrado abaixo.

---

**Isolando $\ddot{\theta}$  na EDO  $I_{zz}\ddot{\theta} + mgr \text{sen}\theta = \tau$, temos**

$$\ddot{\theta}=\frac{1}{I_{zz}}\cdot(\tau-mgr\text{sen}\theta) $$

**Definido $\dot{\theta}=\omega$, temos duas EDO de 1ª ordem que podem ser resolvidas pelas condições inicias fornecidas**

$$\dot{\omega}=\frac{1}{I_{zz}}\cdot(\tau-mgr\text{sen}\theta)$$

$$\dot{\theta}=\omega$$

---

Agora que a equação diferencial já está manipulada de forma apropriada para a solução, podemos calcular o seu resultado. Para verificar se o método está funcionando corretamente, foi simulado e plotado os resultados das soluções das EDO considerando as condições iniciais dadas e $\tau=0$, ou seja, o manipulador deve agir como um pêndulo. Os resultados estão ilustrados abaixo.

![Animação sem controle](https://i.loli.net/2018/07/05/5b3d987d06deb.gif)

Essa imagem mostra os resultados da simulação do manipulador. Como esperado, o comportamento é idêntico ao de um pêndulo simples, movimento harmônico. Isso fica ainda mais claro nos gráficos da posição e velocidade angular abaixo.

![Gráficos](https://i.loli.net/2018/07/05/5b3d99b1ddbd2.png)

Com esses resultados podemos concluir que o método de Range-Kutta de quarta ordem implementado pela biblioteca está funcionando corretamente.

## 2. Projeto do Controlador

O controlador utiliza a técnica de linearização por realimentação, para que o erro decaia de forma exponencial. Para tal, a atuação na dinâmica deve ocorrer através da variável $\tau$. Assim, observando a equaçao diferencial do modelo, podemos implementar o método.

---

**Reorganizado a EDO**

$$\ddot{\theta}=-\frac{1}{I_{zz}}\cdot mgr\text{sen}\theta+\frac{1}{I_{zz}}\cdot\tau $$

**Do método de linearização por realimentação, temos**

$$\ddot{x}=f(x,\dot{x})+b(x,\dot{x})u$$

$$u=b^{-1}(-f(x,\dot{x})+\ddot{x_{d}}-2\lambda \dot{\tilde{x}}-\lambda^2\tilde{x})$$

**Para esse modelo, ficamos com**

$$f(x,\dot{x})=\frac{1}{I_{zz}}\cdot mgr\text{sen}\theta$$

$$b(x,\dot{x})=\frac{1}{I_{zz}}$$

**Portanto, o torque para realizar o controle é**

$$\tau=mgr\text{sen}\theta+I_{zz}(\ddot{\theta_d}-2\lambda \dot{\tilde{\theta}}-\lambda^2\tilde{\theta})$$

---

Assim, só precisamos implementar essa equação no algoritmo de controle para a dada variação no valor de referência $\theta_d$. Após essa implementação, vamos verificar com o comportamento para algumas posições desejadas. A não ser que espefcificado, o valor de $\lambda$ foi tomado como sendo $1$.

### Teste 1 - Posição Angular Constante

O primeiro teste consistiu em verificar a convergência do manipulador para uma posição inicial desejada fixa que correspondeu a $\theta_d=\pi/6$. Nesse caso, a posição pretendida é constante e, portanto, $\dot{\theta_d}=0$ e $\ddot{\theta_d}=0$. Os resultados gráficos e a simulação para esse exemplo estão mostrados abaixo, assumindo condições iniciais nulas $\theta\left(0\right)=0$ e $\dot{\theta}\left(0\right)=0$. 

![5b3dab4b73c51](https://i.loli.net/2018/07/05/5b3dab4b73c51.gif)

A simulação mostra que de fato o algoritmo é capaz de realizar o controle do manipulador de forma satisfatória. Isso fica ainda mais claro observando os gráficos da simulação, onde pode ser notado que a posição angular converge da maneira esperada para o valor esperado.

![5b3daa4cbfdbb](https://i.loli.net/2018/07/05/5b3daa4cbfdbb.png)

Nessa imagem temos vários gráficos interessantes. O primeiro, em vermelho, mostra que o torque inicia em um valor alto e diminui a medida que o manipulador está chegando na posição desejada e se mantém nesse valor. O segundo podemos notar como o erro da velocidade angular faz com que o erro da posição convirja para zero. Por fim, o terceiro e quarto gráficos mostram a dinâmica das variáveis de estados $[\theta,\dot{\theta}]$, seus respectivos erros $[\hat{\theta},\dot{\hat{\theta}}]$, e valores desejados $[\theta_d, \dot{\theta_d}]$. Neles, ficam nítido a convergência da posição angular para o valor desejado.

### Teste 2 - Saturação do Torque

Nesse teste se manteve as condições iniciais do anterior, entretanto, foi definido uma restrição física o motor do manipulador. Essa restrição foi feita através da limitação do torque no intervalo $[-4.8, 4.8]$ , desse modo, exemplificando o limite de torque de um motor real. Então, partindo das mesmas condições inicais anteriores, isto é, $\theta_d=\pi/6$, $\dot{\theta_d}=0$ e $\ddot{\theta_d}=0$, os  resultados gráficos para esse exemplo estão mostrados a seguir.

![5b3db3b38a180](https://i.loli.net/2018/07/05/5b3db3b38a180.gif)

Nesse caso, pode-se notar que a limitação do torque fez com que o controlador fosse incapaz de fazer o manipulador chegar na posição desejada. No caso anterior, o torque máximo aplicado foi de $5 Nm$, e limitando esse valor para apenas $4.8Nm$ já mostrou que o torque é insuficiente. Isso se deu, principalmente, pelo fato do manipulador não ter sido capaz de ultrapasar a posição de 90º que a ação da gravidade no centro de massa causa o maior torque na junta. Assim, podemos perceber que o torque do motor é um dos principais parâmetros para um projeto real de manipulador. 

![5b3db568891ed](https://i.loli.net/2018/07/05/5b3db568891ed.png)

O gráfico do torque revela a limitação em $4.8 Nm$. No segundo gráfico vemos que, diferentemetne do caso anteriro, o erro não converge para $0$. Nos gráficos das variáveis de estados podemos ver o comportamento oscilatório do manipulador causado pela própria dinâmica do movimento. 

## 3. Movimento Oscilatório

Agora para testar se o controlador é capaz de fazer com que o manipulador execute um movimento mais complexo, vamos definir a posição desejada como uma função oscilatória definido pela equação ${\theta}_d = \frac{\pi}{6}\left(1-\text{cos} 2 \pi t \right)$. 

### Teste 3 - Posição Angular Variável

Da primeira e segunda derivada dessa função são retirados os valores de $\dot{\theta}_d$ e $\ddot{\theta}_d$. Então considerando as mesmas condições iniciais utilizadas anteriormente,  $\theta\left(0\right)=0$ e $\dot{\theta}\left(0\right)=0$, os resultados podem ser calculados.

![5b3db93ea8471](https://i.loli.net/2018/07/05/5b3db93ea8471.gif)

Podemos notar que o controlador é capaz de simular a função proposta relativamente rápido. Para agilizar a convergência para a posição desejado poderíamos ajustar o valor de $\lambda$ para outro valor.  Outro ponto interessante que pode ser notado no gráfico é como o período de oscilação é a primiera coisa a ser sincronizada, para que posteriormente a posição seja ajustada aos poucos.

![5b3dbab6a1416](https://i.loli.net/2018/07/05/5b3dbab6a1416.png)

Aqui podemos notar como o comportamento oscilatório da posição desejada reflete um comportamento oscilatório dos gráficos. No gráfico 3, constata-se que a posição está se aproximando cada vez da função desejada ao ponto que, no tempo 10, as curvas já estão praticamente sobrepostas. Também pode-se notar que no gráfico da velocidade angular ocorre uma certa convergência entre a função desejada e  o valor atual nesse intervalo de tempo.

## 4. Adicionando Atrito ao Modelo Simulado

Por fim, iremos adicionar o modelo matemático de atrito seco na equação diferencial da simulação, porém não iremos alterar a equação do controlador. Como explicado em sala de aula, o atrito representa uma incerteza do modelo, de modo que não há garantia se o erro do controlador convergirá para zero. Os testes serão feitos em exemplos que já foram calculados para percebermos as principais diferenças.

### Teste 4 - Refazendo o Teste 1 com Atrito

Todas as condições do Teste 1 foram aproveitadas, no entanto, agora o modelo possui o componente do atrito. Nesse caso, apenas a simulação está mostrada abaixo.

![5b3dc658f2a11](https://i.loli.net/2018/07/05/5b3dc658f2a11.gif)

Pode-se notar que a adição do atrito fez com que o torque calculado pelo controlador fosse insuficiente de ultrapassar a posição de 90º. Assim, fazendo com que o manipulador em equilíbrio na posição inicial, de maneira semelhante ao que aconteceu com o Teste 2, sendo que o torque não está sendo limitado, e sim o controlador que está mandando torque insuficiente. 

Para se ver o desempenho do controlador nesse caso, aumentou-se o valor de $\lambda$ para que o torque necessário para sair dessa posição pudesse ser alcançado. Os resultados com $\lambda=5$ estão mostrados abaixo.

![5b3dc7b2c876d](https://i.loli.net/2018/07/05/5b3dc7b2c876d.gif)

Inicialmente, podemos perceber que o manipulador está se aproximando do valor desejado, entretanto, com o passar do tempo podemos notar que o braço estagna à uma certa distância da posição desejada. Como dito anteriormente, isso é um reflexo de uma parte desconhecida do modelo para o controle e que, consequentemente, não está sendo considerada no controlador. Ou seja, uma das principais desvantagens desse método de controle é a inflexibilidade em relação a incerteza do modelo do manipulador. Isso também faz com que seja um ótimo tipo de controlador em modelos que atuam em ambiente de alto grau de precisão e certeza.

![5b3dc9461c63c](https://i.loli.net/2018/07/05/5b3dc9461c63c.png)

No gráfico do torque podemos notar que ele converge para um valor insuficiente. Nos outros gráficos fica nítido como o ângulo não converge para o valor desejado, apensar de se aproximar consideravelmente.

### Test 5 - Refazendo Teste 3 com Atrito

Mantendo todas as considerações do teste 3, com exceção da adição do atrito e definindo $\lambda=5$. Os resultados estão mostrados abaixo.

![5b3dcab166386](https://i.loli.net/2018/07/05/5b3dcab166386.gif)

Semelhante ao caso anterior, a posição final não nunca é alcançada no tempo desejado. Entretanto, o movimento descrito pela função é alcançado, mesmo que não seja no tempo desejado.

![5b3dcb8cb13ce](https://i.loli.net/2018/07/05/5b3dcb8cb13ce.png)

Nesse caso, vemos que que o erro no final é bem pequeno já que apenas há uma pequena defasagem entre a posição atual e desejada. O gráfico violeta mostra uma geometria estranha causada pela adição do atrito no problema.

# Conclusão

O projeto se mostrou de grande importância no aprofundamento do conhecimento ministrados na disciplina. Desde a análise do comportamento das soluções das equações diferenciais do sistema até a influência que perturbações externas causam no tipo de controle estudado. Também pode-se perceber as vantagens e desvantagens desse tipo de controlador, como a necessidade de um bom modelo dinâmico do sistema.

Por fim, o trabalho também ajudou no desenvolvimento de ferramentas computacionais que podem ser facilmente adaptadas para serem utilizadas em outras disciplinas ou problemas semelhantes.
