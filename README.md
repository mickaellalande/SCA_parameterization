# To implement new Snow Cover Area (SCA) parameterization in LMDZ/ORCHIDEE

Test

The goal is to add the standard deviation in the relationship between Snow Water Equivalent (SWE) and Snow Cover Area (SCA), that is not taken into account so far. The aim would be to reduce the bias over the High Mountain of Asia (HMA) of IPSL-CM6A-LR.

https://docs.google.com/document/d/1gK69TtH3feRFu4q0MjmuouC8xG6Gcth6Qe9cbeY5vIM/edit?usp=sharing

## LMDZOR_v6.1.11
```def
#-H- LMDZOR_v6.1.11  LMDZ and ORCHIDEE model configuration
#-H- LMDZOR_v6.1.11  Configuration corresponding to the coupled modele IPSLCM6.1.11-LR
#-M- LMDZOR_v6.1.11  Josefine.Ghattas@ipsl.fr
#-C- LMDZOR_v6.1.11  IOIPSL/tags/v2_2_4/src             HEAD 8    IOIPSL/src   modeles
#-C- LMDZOR_v6.1.11  tags/ORCHIDEE_2_0/ORCHIDEE         6592 14   ORCHIDEE     modeles
#-C- LMDZOR_v6.1.11  LMDZ6/branches/IPSLCM6.0.15        3643 11   LMDZ         modeles
#-C- LMDZOR_v6.1.11  XIOS/branchs/xios-2.5              1903 12   XIOS         modeles
#-C- LMDZOR_v6.1.11  CONFIG/UNIFORM/v6/LMDZOR_v6        4914 8    LMDZOR_v6    config
#-C- LMDZOR_v6.1.11  tags/ORCHIDEE_2_0/ORCHIDEE_OL      6592 14   ORCHIDEE_OL  config
#-C- LMDZOR_v6.1.11  trunk/libIGCM                      1534 10   libIGCM      .  
```

### LMDZ version
```bash
Working Copy Root Path: /gpfsdswork/projects/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/LMDZ
URL: http://svn.lmd.jussieu.fr/LMDZ/LMDZ6/branches/IPSLCM6.0.15
Relative URL: ^/LMDZ6/branches/IPSLCM6.0.15
Repository Root: http://svn.lmd.jussieu.fr/LMDZ
Repository UUID: e51f81be-29bc-408f-98e3-ee85b5628ff9
Revision: 3616
Node Kind: directory
Schedule: normal
Last Changed Author: oboucher
Last Changed Rev: 3616
Last Changed Date: 2019-12-16 14:10:03 +0100 (Mon, 16 Dec 2019)
```

### ORCHIDEE version
```bash
Working Copy Root Path: /gpfsdswork/projects/rech/goe/ufz23bm/SCA_parameterization/modipsl/modeles/ORCHIDEE
URL: svn://forge.ipsl.jussieu.fr/orchidee/tags/ORCHIDEE_2_0/ORCHIDEE
Relative URL: ^/tags/ORCHIDEE_2_0/ORCHIDEE
Repository Root: svn://forge.ipsl.jussieu.fr/orchidee
Repository UUID: f489ceea-5127-0410-b15c-c4a6149ed9a7
Revision: 6592
Node Kind: directory
Schedule: normal
Last Changed Author: josefine.ghattas
Last Changed Rev: 6592
Last Changed Date: 2020-03-04 10:34:45 +0100 (Wed, 04 Mar 2020)
```

The original code is left on the `master` branch and the developments will be done on the `test` branch. Objectives:
1. Get the std not averaged in `modipsl/modeles/LMDZ/libf/phylmd/grid_noro_m.F90`
2. Use it in `modipsl/modeles/ORCHIDEE/src_sechiba/condveg.f90` for changing the snow cover fraction relationship

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
