# Numerical Solution of an Ordinary Differential Equation in ANSI C

This project implements several numerical methods for solving an initial value problem (IVP) for an ordinary differential equation.

The solved equation is:

dy/dx = xy

with the initial condition:

y(0) = 1

The analytical solution is:

y = exp(x²/2)

## Implemented Methods

- Euler Method
- Improved Euler Method (Heun's Method)
- Classical Runge-Kutta Method of 4th order (RK4)

## Features

- Comparison between numerical and analytical solutions
- Calculation of absolute errors
- Error analysis for different step sizes
- Data export for plotting
- Visualization using Gnuplot

## Technologies

- ANSI C
- GCC Compiler
- Gnuplot

## Purpose

The goal of this project is to study the accuracy and behavior of different numerical methods for solving differential equations.
