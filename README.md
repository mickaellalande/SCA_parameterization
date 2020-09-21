# To implement new Snow Cover Area (SCA) parameterization in LMDZ/ORCHIDEE

The goal is to add the standard deviation in the relationship between Snow Water Equivalent (SWE) and Snow Cover Area (SCA), that is not taken into account so far. The aim would be to reduce the bias over the High Mountain of Asia (HMA) of IPSL-CM6A-LR.

https://docs.google.com/document/d/1gK69TtH3feRFu4q0MjmuouC8xG6Gcth6Qe9cbeY5vIM/edit?usp=sharing

## Bibliography

### Gerhard implementation
A new snow cover fraction parametrization for the ECHAM4 GCM ([Roesch et al., 2001](https://link.springer.com/article/10.1007/s003820100153))
- https://github.com/mickaellalande/PhD/blob/master/local/SCE_SWE_parametization/Roesch2001.ipynb

### What is in Orchidée (https://orchidas.lsce.ipsl.fr/dev/albedo/)
An observation‐based formulation of snow cover fraction and its evaluation over large North American river basins ([Niu and Yang, 2007](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2007JD008674))

-  https://github.com/mickaellalande/PhD/blob/master/local/SCE_SWE_parametization/Niu2007.ipynb  
- https://github.com/mickaellalande/PhD/blob/master/local/SCE_SWE_parametization/Niu2007-std.ipynb

### New version with arccos
A new fractional snow‐covered area parameterization for the Community Land Model and its effect on the surface energy balance ([Swenson and Lawrence, 2012](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2012JD018178))  

- https://github.com/mickaellalande/PhD/blob/master/local/SCE_SWE_parametization/Swenson2012.ipynb
