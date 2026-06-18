# Nonlinear Programming Assignment

Final assignment for the course **Optimization for Systems and Control (SC42056)** at TU Delft.

This project investigates freeway congestion control using ramp metering optimization based on the discrete-time macroscopic traffic flow model **METANET**. Different nonlinear optimization approaches are implemented and compared in MATLAB to minimize the Total Time Spent (TTS) of vehicles on a freeway segment.

The assignment explores:

- Nonlinear constrained optimization
- Sequential Quadratic Programming (SQP)
- Interior Point (IP) optimization
- Genetic Algorithms (GA)
- Dynamic traffic flow modelling
- Ramp metering control strategies


## Project Overview

The freeway system is modelled using the METANET traffic model, where the objective is to optimize the ramp metering rate \( r(k) \) in order to reduce congestion and minimize the Total Time Spent (TTS) by drivers.

Three optimization approaches are investigated throughout the assignment and named as:

- **BigMac**
- **RudiMINtal**
- **Rudi**

The repository contains MATLAB implementations of the optimization algorithms together with the [final report](NLP_report.pdf) containing derivations, simulations, plots, and analysis results.

## Repository Structure

```text
Nonlinear-Programming-Ramp-Metering/
│
├── Code/
│   ├── bigmac.m
│   ├── rudi.m
│   └── rudimintal_and_discrete.m
│
├── NLP_report.pdf
│
└── README.md
```
---

## Collaborators

- Melis Orhun
- Sven Rutgers
