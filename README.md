# astrocytes_population_sim

## General Info

This is the code of the paper "Engineering Molecular Communications of Astrocytes for Logic Gate Operations". You may find the paper at https://arxiv.org/abs/2007.06646.

We use the SSA Guillespie method for solving the individual astrocyte calcium signalling ODE over a population. Please refer to these papers for more info

- M.T. Barros, S. Balasubramaniam, B. Jennings. Comparative End-to-end Analysis of Ca2+ Signaling-based Molecular Communication in Biological Tissues. IEEE Transactions on Communications, v. 63, p 5128-5142. 2015.
- M. Lavrentovich and S. Hemkin, “A mathematical model of spontaneous calcium(II) oscillations in astrocytes,” J. Theor. Biol., vol. 251, pp. 553– 560, 2008.
- D. T. Gillespie, "Exact stochastic simulation of coupled chemical reactions" Journal of Physical Chemistry, vol. 81, no. 25, pp. 2340-2361, 1977.

## How to use it

You may run the code with initial configuration
- python astrocytes_population_sim.py
  
Alternatevely, the code supports the pypy interpreter, therefore use
- pypy astrocytes_population_sim.py

## Notes

- The code available will not automatically reproduce the content of the paper. It will give the code of the simulator used, however, you may need to understand each graph in order to adjust the simulator that can give the same results of the paper. 
- This code will only work for the theoretical results only!!
