Author information

-----------------------------------------------

Yuan-Chao Hu (ychu0213@gmail.com)

SSLAB

--------------------------------------------------------



This module provide LJ type Dzugutov potential for MD simulations,
the form of the potential is :

$$
V = V_1 + V_2
$$

$$
V_1 = A(r^{-m} - B) \exp(\frac{c}{r-a}) \qquad (r \lt a)
$$

$$
V_2 = B \exp(\frac{d}{r - b}) \qquad (r \lt b)
$$

The parameters have the following values:




The usage from lammps is as following"

----

pair_style       dzugutov 1.94

pair coeff       1 1 1.94

----

The last number for each line is the cutoff distance.
The other parameters are included in the module itself.


Benchmarking of analytical function to lammps simulations:

![compare](./test/twoatoms/compare_twoatoms.png)