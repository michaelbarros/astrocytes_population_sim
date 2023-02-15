# astrocytes_population_sim

## General Info

This is the code of the paper "Engineering Molecular Communications of Astrocytes for Logic Gate Operations". You may find the paper at https://arxiv.org/abs/2007.06646.

We use the SSA Guillespie method for solving the individual astrocyte calcium signalling ODE over a population. Please refer to the references below for more detailed information.

## How to use it

You may run the code with initial configuration
- python astrocytes_population_sim.py
  
Alternatevely, the code supports the pypy interpreter, therefore use
- pypy astrocytes_population_sim.py

# Minimum Requirements

- Python 2.7 and numpy

# Code

We used Guillespie exact stochastic method alongside with the mininal model for calcium wave oscillation. If the user wants further knowledge about the simulation.

Parameter   | Description
---------   | ------
destination | The destination cell
conc        | Concentration of Calcium on Receiver
n, m        | 3D Tissue Size
x           | Length of the cell
txconc      | Concentration of Calcium on Transmitter
deltaamp    | Amplification of Calcium

## Other References

1. [M. T. Barros, S. Balasubramaniam, and B. Jennings, Comparative End-to-end Analysis of Ca2+ Signaling-based Molecular Communication in Biological Tissues. IEEE Transactions on Communications, v. 63, p 5128-5142. 2015.](https://sites.google.com/site/michaeltaob/07289389.pdf?attredirects=0&d=1)

2. [M. T. Barros, S. Balasubramaniam, B. Jennings, and Y. Koucheryavy, Transmission protocols for calcium signaling based molecular communications in deformable cellular tissues, Manuscript sumbitted for journal publication, 2014.](https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxtaWNoYWVsdGFvYnxneDo0NzMyNGNkOTFiOGZjZjll)

3. [M. T. Barros, S. Balasubramaniam, and B. Jennings, Integrating information metrics and molecular communications for inferencing and detecting cellular tissue deformation,  IEEE Transactions on Nanobioscience, v. 13, p. 278-288, 2014.](https://docs.google.com/viewer?a=v&pid=sites&srcid=ZGVmYXVsdGRvbWFpbnxtaWNoYWVsdGFvYnxneDo3YzgxYWZlNTFmZDA4ZjQz)

4. [A. Goldbeter, G. Dupont, and M. J. Berridge,  Minimal model for signalinduced Ca2+ oscillations and for their frequency encoding through protein phosphorylation, Proceedings of the National Academy of Science USA, vol. 87, pp. 1461-1465, 1990.](https://www.ncbi.nlm.nih.gov/pubmed/2304911)

5.  [D. T. Gillespie, Exact stochastic simulation of coupled chemical reactions, Journal of Physical Chemistry, vol. 81, no. 25, pp. 2340-2361, 1977.](http://pubs.acs.org/doi/abs/10.1021/j100540a008?journalCode=jpchax)

## Notes

- The code available will not automatically reproduce the content of the paper. It will give the code of the simulator used, however, you may need to understand each graph in order to adjust the simulator that can give the same results of the paper. 
- This code will only work for the theoretical results only!!
