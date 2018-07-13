# Robótica - Projetos

Pequena coleção de projetos desenvolvidos ao longo do Mestrado em Engenharia Mecatrônica da UFRN. Os problemas vão desde a solução da cinemática de alguns manipuladores usando `Sympy`, até a simulação de um controlador para a dinâmica de uma manipulador simples.

### Requisitos

* Sympy (Cálculo da cinemática e dinâmica de forma algébrica)

* Matplolib (Visualização da simulação)

* Scipy (Resolução das equações dinâmicas)

* Dill (Para salvar / recuperar resultados já calculados.)

## Sumário

1. [Cinemática](#cinemática)

2. [Dinâmica](#dinâmica)

3. [Controle](#controle)

### [Cinemática](https://github.com/Fernandohf/Robotica-Projetos/tree/master/1-Cinem%C3%A1tica)

Essa pasta contém Notebooks e Scripts que exemplicam o raciocínio da utilização do [Sympy](http://www.sympy.org/pt/index.html) na resolução da cinemática de alguns manipuladores clássicos. Os exemplos dados abordam desde a definição das variáveis, até a solução da cinemática dos pontos de forma numérica.

### [Dinâmica](https://github.com/Fernandohf/Robotica-Projetos/tree/master/2-Din%C3%A2mica)

Como nos aquivos da cinemática, porém esses Notebooks exemplicam o desenvolvimento das equações de movimento dos manipuladores usando esse o método de Lagrange. Nesse caso, a biblioteca utilizada para esses cálculos ainda é o [Sympy](http://www.sympy.org/pt/index.html), mais especficicamente o módulo de Fisica ([Physics Module](http://docs.sympy.org/latest/modules/physics/index.html)).

### [Controle](https://github.com/Fernandohf/Robotica-Projetos/tree/master/3-Controle)

Dada a dinâmica de um manipulador, isto é, sua equação de movimento, o seu controle é implementado através do método de linearização por realimentação. Os resultados são simulados e exibidos em gráficos da [Matplotlib](https://matplotlib.org/).
