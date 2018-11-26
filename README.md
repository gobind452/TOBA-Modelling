# TOBA-Modelling

This repository contains code used for data analysis and modelling an optical lever attached to a torsion pendulum. This model is used for predicting noise and signals obtained from TOBA. 

Read more about TOBA here - https://journals.aps.org/prd/abstract/10.1103/PhysRevD.95.082004

## Schematics
+ Lagrangian.py - Model the dynamics of an asymmetric torsion pendulum
+ calculateQ.py - Analyze the data obtained to calculate Q(Quality Factor)
+ getReflectionY - Get the signal obtained from the optical lever

Rest are helper functions.

## Issues 
+ Lots of hardcoding
+ Lots of tunable parameters
