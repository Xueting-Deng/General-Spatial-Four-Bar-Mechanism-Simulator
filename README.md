# Spatial Four-Bar Mechanism - Simulation and Database Generation

<p align="center">
<img width="300" height="200" alt="image" src="https://github.com/user-attachments/assets/2012115c-b3ce-4392-9c0e-98c85b7c8e7f" />
<img width="278" height="278" alt="Picture1" src="https://github.com/user-attachments/assets/5c614f2e-7a28-4606-b458-8c8e5f60f252" />
</p>

## Introduction
This repository contains the source code for a general spatial four-bar mechanism simulator and the corresponding database generation pipeline developed for our paper (the paper link will be added once published). For details about the representing mechanism, simulation pipeline, and database generation, please refer to the paper as well as Xueting's phd dissertation (links will be attached here after publication).

The current version of the simulator explicitly encodes all 26 topologies of 1-DOF spatial four-bar mechanisms in the `mec_topo.py` file. While the framework is extensible, this released version may still contain unresolved issues when applied to other mechanism topologies, including spatial six-bar mechanisms.

For access to the newer and more generalized version of the simulator, please contact our lab advisor: Prof. Anurag Purwar (anurag.purwar@stonybrook.edu).

## Required packages
- Python (tested on 3.7.1 / 3.12.4)
- Numpy (tested on 1.21.5 / 2.1.2)
- Scipy (tested on 1.7.3 / 1.13.0)
- Sympy (tested on 1.10.1 / 1.13.3)
- Matplotlib (tested on 3.5.3)
- Vpython (tested on 7.6.4)

## Table of Files
├── Simulator - Jupyter Notebook.ipynb #example code of simulating a mechanism  
├── Simulator.py #example code of simulating a mechanism  
├── Simulator - heart.py #example code of simulating a mechanism that can generate a "heart" shape curve  
├── Simulator_Functions.py #basic functions for simulator  
├── Simulator_Generate_Database.py #example script for generating a database 
├── mec_topo.py #records the basic 26 1-DOF spatial four-bar mechanisms' topology to help the simulator and database generation 

![](https://github.com/Xueting-Deng/General-Spatial-Four-Bar-Mechanism-Simulator/blob/main/heart.gif)

