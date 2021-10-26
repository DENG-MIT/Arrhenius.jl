# Arrhenius

We are in an early-development. Expect some adventures and rough edges.

## Installation

> pkg> add https://github.com/DENG-MIT/Arrhenius.jl


## Publication

+ [Arrhenius.jl: A Differentiable Combustion Simulation Package](https://www.researchgate.net/publication/350573212_Arrheniusjl_A_Differentiable_Combustion_Simulation_Package): overview of Arrhenius.jl and applications in deep mechanism reduction, uncertainty quantification, mechanism tuning and model discovery. [Slides in NCM21](https://www.slideshare.net/WeiqiJi/arrheniusjl-a-differentiable-combustion-simulation-package-248457895), [Vedio for NCM21](https://www.youtube.com/watch?v=X1mwpW78NvA).
+ [Machine Learning Approaches to Learn HyChem Models](https://www.researchgate.net/publication/350890609_Machine_Learning_Approaches_to_Learn_HyChem_Models): demonstrate 1000 times faster than genetic algorithms using commercial software for optimizing complex kinetic models.
+ [Neural Differential Equations for Inverse Modeling in Model Combustors](https://www.researchgate.net/publication/351223124_Neural_Differential_Equations_for_Inverse_Modeling_in_Model_Combustors)


## Applications

+ **Sensitivity analysis for auto-ignition** | [repo](https://github.com/DENG-MIT/ArrheniusActiveSubspace) | Features: auto-differentiation, multi-threading, sensitivity to all of three Arrhenius params A, b and Ea, active subspace based uncertainty quantification
+ **Sensitivity analysis for one-dimensional flames** | [repo](https://github.com/DENG-MIT/Arrhenius_Flame_1D) | Features: auto-differentiation, multi-threading, sensitivity to all of three Arrhenius params A, b and Ea.
+ **Automonous learn kinetic mechanism using neural network** | [repo](https://github.com/DENG-MIT/CRNN_HyChem) | Features: Chemical Reaction Neural Network (CRNN), Neural Ordinary Differential Equations.
+ **Deep Reduction** | [repo](https://github.com/DENG-MIT/DeepReduction) | Features: Two-stages mechanism reduction with deep learning.

**Examples**

> Note that some of the examples are in development and you can have early access by contacting [Weiqi Ji](mailto:weiqiji@mit.edu)
  + [Pyrolysis of JP10](./example/pyrolysis/pyrolysis.ipynb)
  + [Perfect Stirred Reactor](./example/perfect_stirred_reactor)
  + [Auto-ignition](https://github.com/DENG-MIT/NN-Ignition):
  + [Compute Jacobian using AD](https://gist.github.com/jiweiqi/21b8d149bd95b97d9ae948ab92e446df)
