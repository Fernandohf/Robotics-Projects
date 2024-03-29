# Robotics

Collections of Robotics projects. The notebooks covers from kinematics of simple manipulators, to control simulations using **Python**.

## Robots

### Scara

![Scara robot motion animation](https://github.com/Fernandohf/Robotics-Projects/blob/master/media/scara.gif?raw=true)

Scara is a robotic manipulator with 2 degrees of freedom (DoF) commonly used for pick and place tasks.

### Anthropomorphic

![Anthropomorphic arm motion animation](https://github.com/Fernandohf/Robotics-Projects/raw/master/media/msedge_xbxfGl3dEd.gif)

Anthropomorphic robots are inspired in humam arms and can serve multiple purpose tasks. The version studied here has 3 DoF.

## Requirements

- Sympy (Kinematics e Dynamics)
- Matplolib (Visualization)
- Scipy (Dynamics Equations Solver)
- Dill (Save / Load Previous Results)

## Summary

1. [Kinematics](#kinematics)
2. [Dynamics](#dynamics)
3. [Control](#control)

## Kinematics

Kinematics studies the relationship between the joints angles and how the produce desire movement on multiple degree of freedom linkages. This folder contains a collection of notebooks and scripts explaining how to use [Sympy](http://www.sympy.org/pt/index.html) to solve kinematics of 2 robots: [Scara](#scara) and [Anthropomorphic Arm](#anthropomorphic).

## Dynamics

These Notebooks explain how to develop movement equations from robotics arm, using Lagrange's method from [Sympy - Physics Module](http://docs.sympy.org/latest/modules/physics/index.html).

## Simulation

Here the models is simulated using `matplotlib` animations. Additionally, the resulting dynamics equations are reused using `pickle`.

![Pendulum limulation animation](https://github.com/Fernandohf/Robotics-Projects/blob/2353144a9388bebd67827767bd128dd93b46ffd0/media/animation-output.gif)

## Control

Given previous results, it is now implemented linearization by state feedback. The results are animated and visualized using [Matplotlib](https://matplotlib.org/).

### Contents

- Linearization by State Feedback
