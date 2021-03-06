!OCL SCALAR
SUBROUTINE RRTM_KGB4_00

!     Originally by Eli J. Mlawer, Atmospheric & Environmental Research.
!     BAND 4:  630-700 cm-1 (low - H2O,CO2; high - O3,CO2)
!     Reformatted for F90 by JJMorcrette, ECMWF
!     Reformatted for NEC by H.Gall�e   , LGGE  (splitting)

!     ------------------------------------------------------------------

#include "tsmbkind.h"

USE YOERRTO4 , ONLY : KAO     ,KBO     ,SELFREFO   ,FRACREFAO  ,FRACREFBO
USE YOERRTA4 , ONLY : STRRAT1   ,STRRAT2

!     ------------------------------------------------------------------


IMPLICIT NONE
STRRAT1 = 850.577_JPRB
STRRAT2 = 35.7416_JPRB

!     ------------------------------------------------------------------


!     The array SELFREFO contains the coefficient of the water vapor
!     self-continuum (including the energy term).  The first index
!     refers to temperature in 7.2 degree increments.  For instance,
!     JT = 1 refers to a temperature of 245.6, JT = 2 refers to 252.8,
!     etc.  The second index runs over the g-channel (1 to 16).

SELFREFO( :, 1) = (/&
&2.62628E-01_JPRB, 2.29008E-01_JPRB, 1.99692E-01_JPRB, 1.74129E-01_JPRB, 1.51838E-01_JPRB,&
&1.32400E-01_JPRB, 1.15451E-01_JPRB, 1.00672E-01_JPRB, 8.77845E-02_JPRB, 7.65469E-02_JPRB/)
SELFREFO( :, 2) = (/&
&2.45051E-01_JPRB, 2.12961E-01_JPRB, 1.85073E-01_JPRB, 1.60837E-01_JPRB, 1.39775E-01_JPRB,&
&1.21471E-01_JPRB, 1.05564E-01_JPRB, 9.17397E-02_JPRB, 7.97260E-02_JPRB, 6.92856E-02_JPRB/)
SELFREFO( :, 3) = (/&
&2.42194E-01_JPRB, 2.09976E-01_JPRB, 1.82044E-01_JPRB, 1.57827E-01_JPRB, 1.36832E-01_JPRB,&
&1.18630E-01_JPRB, 1.02849E-01_JPRB, 8.91673E-02_JPRB, 7.73057E-02_JPRB, 6.70221E-02_JPRB/)
SELFREFO( :, 4) = (/&
&2.44485E-01_JPRB, 2.11926E-01_JPRB, 1.83702E-01_JPRB, 1.59237E-01_JPRB, 1.38030E-01_JPRB,&
&1.19648E-01_JPRB, 1.03714E-01_JPRB, 8.99014E-02_JPRB, 7.79286E-02_JPRB, 6.75503E-02_JPRB/)
SELFREFO( :, 5) = (/&
&2.43120E-01_JPRB, 2.10743E-01_JPRB, 1.82679E-01_JPRB, 1.58351E-01_JPRB, 1.37263E-01_JPRB,&
&1.18984E-01_JPRB, 1.03139E-01_JPRB, 8.94038E-02_JPRB, 7.74978E-02_JPRB, 6.71774E-02_JPRB/)
SELFREFO( :, 6) = (/&
&2.40558E-01_JPRB, 2.08922E-01_JPRB, 1.81446E-01_JPRB, 1.57583E-01_JPRB, 1.36859E-01_JPRB,&
&1.18860E-01_JPRB, 1.03229E-01_JPRB, 8.96529E-02_JPRB, 7.78624E-02_JPRB, 6.76225E-02_JPRB/)
SELFREFO( :, 7) = (/&
&2.42496E-01_JPRB, 2.10386E-01_JPRB, 1.82528E-01_JPRB, 1.58359E-01_JPRB, 1.37390E-01_JPRB,&
&1.19198E-01_JPRB, 1.03415E-01_JPRB, 8.97211E-02_JPRB, 7.78409E-02_JPRB, 6.75337E-02_JPRB/)
SELFREFO( :, 8) = (/&
&2.39781E-01_JPRB, 2.08227E-01_JPRB, 1.80825E-01_JPRB, 1.57029E-01_JPRB, 1.36365E-01_JPRB,&
&1.18419E-01_JPRB, 1.02836E-01_JPRB, 8.93030E-02_JPRB, 7.75510E-02_JPRB, 6.73456E-02_JPRB/)
SELFREFO( :, 9) = (/&
&2.38707E-01_JPRB, 2.07058E-01_JPRB, 1.79605E-01_JPRB, 1.55792E-01_JPRB, 1.35136E-01_JPRB,&
&1.17219E-01_JPRB, 1.01677E-01_JPRB, 8.81962E-02_JPRB, 7.65026E-02_JPRB, 6.63594E-02_JPRB/)
SELFREFO( :,10) = (/&
&2.29942E-01_JPRB, 2.00668E-01_JPRB, 1.75121E-01_JPRB, 1.52826E-01_JPRB, 1.33370E-01_JPRB,&
&1.16390E-01_JPRB, 1.01572E-01_JPRB, 8.86410E-02_JPRB, 7.73560E-02_JPRB, 6.75077E-02_JPRB/)
SELFREFO( :,11) = (/&
&2.39870E-01_JPRB, 2.08120E-01_JPRB, 1.80573E-01_JPRB, 1.56671E-01_JPRB, 1.35934E-01_JPRB,&
&1.17941E-01_JPRB, 1.02330E-01_JPRB, 8.87854E-02_JPRB, 7.70335E-02_JPRB, 6.68371E-02_JPRB/)
SELFREFO( :,12) = (/&
&2.40196E-01_JPRB, 2.08400E-01_JPRB, 1.80812E-01_JPRB, 1.56877E-01_JPRB, 1.36110E-01_JPRB,&
&1.18092E-01_JPRB, 1.02460E-01_JPRB, 8.88962E-02_JPRB, 7.71284E-02_JPRB, 6.69184E-02_JPRB/)
SELFREFO( :,13) = (/&
&2.40426E-01_JPRB, 2.08603E-01_JPRB, 1.80991E-01_JPRB, 1.57035E-01_JPRB, 1.36249E-01_JPRB,&
&1.18214E-01_JPRB, 1.02567E-01_JPRB, 8.89909E-02_JPRB, 7.72117E-02_JPRB, 6.69917E-02_JPRB/)
SELFREFO( :,14) = (/&
&2.40590E-01_JPRB, 2.08742E-01_JPRB, 1.81110E-01_JPRB, 1.57135E-01_JPRB, 1.36334E-01_JPRB,&
&1.18287E-01_JPRB, 1.02628E-01_JPRB, 8.90428E-02_JPRB, 7.72556E-02_JPRB, 6.70288E-02_JPRB/)
SELFREFO( :,15) = (/&
&2.40634E-01_JPRB, 2.08779E-01_JPRB, 1.81141E-01_JPRB, 1.57162E-01_JPRB, 1.36357E-01_JPRB,&
&1.18306E-01_JPRB, 1.02645E-01_JPRB, 8.90565E-02_JPRB, 7.72673E-02_JPRB, 6.70387E-02_JPRB/)
SELFREFO( :,16) = (/&
&2.40652E-01_JPRB, 2.08793E-01_JPRB, 1.81151E-01_JPRB, 1.57169E-01_JPRB, 1.36362E-01_JPRB,&
&1.18309E-01_JPRB, 1.02647E-01_JPRB, 8.90576E-02_JPRB, 7.72675E-02_JPRB, 6.70383E-02_JPRB/)

FRACREFAO( :, 1) = (/&
!     From P = 
    &0.15579100_JPRB,0.14918099_JPRB,0.14113800_JPRB,0.13127001_JPRB,&
    &0.11796300_JPRB,0.10174300_JPRB,0.08282370_JPRB,0.06238150_JPRB,&
    &0.04213440_JPRB,0.00458968_JPRB,0.00377949_JPRB,0.00298736_JPRB,&
    &0.00220743_JPRB,0.00140644_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 2) = (/&
    &0.15292799_JPRB,0.15004000_JPRB,0.14211500_JPRB,0.13176700_JPRB,&
    &0.11821100_JPRB,0.10186300_JPRB,0.08288040_JPRB,0.06241390_JPRB,&
    &0.04220720_JPRB,0.00459006_JPRB,0.00377919_JPRB,0.00298743_JPRB,&
    &0.00220743_JPRB,0.00140644_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 3) = (/&
    &0.14386199_JPRB,0.15125300_JPRB,0.14650001_JPRB,0.13377000_JPRB,&
    &0.11895900_JPRB,0.10229400_JPRB,0.08312110_JPRB,0.06239520_JPRB,&
    &0.04225560_JPRB,0.00459428_JPRB,0.00378865_JPRB,0.00298860_JPRB,&
    &0.00220743_JPRB,0.00140644_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 4) = (/&
    &0.14359100_JPRB,0.14561599_JPRB,0.14479300_JPRB,0.13740200_JPRB,&
    &0.12150100_JPRB,0.10315400_JPRB,0.08355480_JPRB,0.06247240_JPRB,&
    &0.04230980_JPRB,0.00459916_JPRB,0.00378373_JPRB,0.00300063_JPRB,&
    &0.00221111_JPRB,0.00140644_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 5) = (/&
    &0.14337599_JPRB,0.14451601_JPRB,0.14238000_JPRB,0.13520500_JPRB,&
    &0.12354200_JPRB,0.10581200_JPRB,0.08451810_JPRB,0.06262440_JPRB,&
    &0.04239590_JPRB,0.00460297_JPRB,0.00378701_JPRB,0.00300466_JPRB,&
    &0.00221899_JPRB,0.00141020_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 6) = (/&
    &0.14322001_JPRB,0.14397401_JPRB,0.14117201_JPRB,0.13401900_JPRB,&
    &0.12255500_JPRB,0.10774100_JPRB,0.08617650_JPRB,0.06296420_JPRB,&
    &0.04249590_JPRB,0.00463406_JPRB,0.00378241_JPRB,0.00302037_JPRB,&
    &0.00221583_JPRB,0.00141103_JPRB,0.00053814_JPRB,0.00007991_JPRB/)
FRACREFAO( :, 7) = (/&
    &0.14309500_JPRB,0.14364301_JPRB,0.14043900_JPRB,0.13348100_JPRB,&
    &0.12211600_JPRB,0.10684700_JPRB,0.08820590_JPRB,0.06374610_JPRB,&
    &0.04264730_JPRB,0.00464231_JPRB,0.00384022_JPRB,0.00303427_JPRB,&
    &0.00221825_JPRB,0.00140943_JPRB,0.00055564_JPRB,0.00007991_JPRB/)
FRACREFAO( :, 8) = (/&
    &0.15579100_JPRB,0.14918099_JPRB,0.14113800_JPRB,0.13127001_JPRB,&
    &0.11796300_JPRB,0.10174300_JPRB,0.08282370_JPRB,0.06238150_JPRB,&
    &0.04213440_JPRB,0.00458968_JPRB,0.00377949_JPRB,0.00298736_JPRB,&
    &0.00220743_JPRB,0.00140644_JPRB,0.00053024_JPRB,0.00007459_JPRB/)
FRACREFAO( :, 9) = (/&
    &0.15937001_JPRB,0.15159500_JPRB,0.14242800_JPRB,0.13078900_JPRB,&
    &0.11671300_JPRB,0.10035700_JPRB,0.08143450_JPRB,0.06093850_JPRB,&
    &0.04105320_JPRB,0.00446233_JPRB,0.00369844_JPRB,0.00293784_JPRB,&
    &0.00216425_JPRB,0.00143403_JPRB,0.00054571_JPRB,0.00007991_JPRB/)

FRACREFBO( :, 1) = (/&
!     From P = 1.17 mb.
    &0.15558299_JPRB,0.14930600_JPRB,0.14104301_JPRB,0.13124099_JPRB,&
    &0.11792900_JPRB,0.10159200_JPRB,0.08314130_JPRB,0.06240450_JPRB,&
    &0.04217020_JPRB,0.00459313_JPRB,0.00379798_JPRB,0.00299835_JPRB,&
    &0.00218950_JPRB,0.00140615_JPRB,0.00053010_JPRB,0.00007457_JPRB/)
FRACREFBO( :, 2) = (/&
    &0.15592700_JPRB,0.14918999_JPRB,0.14095700_JPRB,0.13115700_JPRB,&
    &0.11788900_JPRB,0.10158000_JPRB,0.08313780_JPRB,0.06240240_JPRB,&
    &0.04217000_JPRB,0.00459313_JPRB,0.00379798_JPRB,0.00299835_JPRB,&
    &0.00218950_JPRB,0.00140615_JPRB,0.00053010_JPRB,0.00007457_JPRB/)
FRACREFBO( :, 3) = (/&
    &0.15949000_JPRB,0.15014900_JPRB,0.14162201_JPRB,0.13080800_JPRB,&
    &0.11713500_JPRB,0.10057100_JPRB,0.08170080_JPRB,0.06128110_JPRB,&
    &0.04165600_JPRB,0.00459202_JPRB,0.00379835_JPRB,0.00299717_JPRB,&
    &0.00218958_JPRB,0.00140616_JPRB,0.00053010_JPRB,0.00007457_JPRB/)
FRACREFBO( :, 4) = (/&
    &0.15967900_JPRB,0.15038200_JPRB,0.14196999_JPRB,0.13074800_JPRB,&
    &0.11701700_JPRB,0.10053000_JPRB,0.08160790_JPRB,0.06122690_JPRB,&
    &0.04128310_JPRB,0.00456598_JPRB,0.00379486_JPRB,0.00299457_JPRB,&
    &0.00219016_JPRB,0.00140619_JPRB,0.00053011_JPRB,0.00007456_JPRB/)
FRACREFBO( :, 5) = (/&
    &0.15989800_JPRB,0.15057300_JPRB,0.14207700_JPRB,0.13068600_JPRB,&
    &0.11682900_JPRB,0.10053900_JPRB,0.08163610_JPRB,0.06121870_JPRB,&
    &0.04121690_JPRB,0.00449061_JPRB,0.00371235_JPRB,0.00294207_JPRB,&
    &0.00217778_JPRB,0.00139877_JPRB,0.00053011_JPRB,0.00007455_JPRB/)
FRACREFBO( :, 6) = (/&
    &0.15950100_JPRB,0.15112500_JPRB,0.14199100_JPRB,0.13071300_JPRB,&
    &0.11680800_JPRB,0.10054600_JPRB,0.08179050_JPRB,0.06120910_JPRB,&
    &0.04126050_JPRB,0.00444324_JPRB,0.00366843_JPRB,0.00289369_JPRB,&
    &0.00211550_JPRB,0.00134746_JPRB,0.00050874_JPRB,0.00007863_JPRB/)

!     -----------------------------------------------------------------
RETURN
END SUBROUTINE RRTM_KGB4_00
