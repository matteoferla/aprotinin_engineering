# aprotinin_engineering
Improve the interface between aprotinin and Dengue protease in structure PDB:3U1J.

This interface is partially fortuitous as aprotinin, also known as bovine pancreatic trypsin inhibitor, is a broad spectrum serine protease inhibitor. Therefore the is likely room for improvement.
The reason why this was requested was told to me but I forgot.
That the active site is occluded is not a problem I am told.

## Plan

* Minimize
* Point mutation scan
* FastDesign without cysteines
* FastDesign interface

## Minimization

Standard, coordinate restrained

> Code: [01-minimise.py](01-minimise.py)

## Cysteines

The protein has 3 disulfide bonds (C5-C55, C30-C51, C14-C38). There are a few designs were cysteine to alanine variants are made, such as PDB:2ZJX. Cysteine to alanine is conservative, but might not be the best.

> Code: [02-remove_cystine.py](02-remove_cystine.py)

Rosetta Relax with design found these pairs were good:

* C5V C55A
* C14A C38V
* C30A C51T

Together the predicted abs ∆G goes from -745 to -764 kcal/mol.	

## Point mutations

> Code: [03-scan.py](03-scan.py)

> NB. All designs were on the relaxed wild type, not a C5V C14A C30A C38V C51T C55A variant


| mutations | complex_ddG |   interface_ddG | 
|:----------|------------:|----------------:|
| I18L      |         1.4 |            -2.6 | 
| I19W      |         1.6 |            -2.1 | 
| G36G*     |        -0.5 |            -2   |  
| P13E      |         1.7 |            -1.5 |
| R17M      |         1.1 |            -1.4 | 
| R17L      |        -0.2 |            -1.3 | 
| R17Q      |         0.9 |            -0.7 | 
| P13N      |         1.3 |            -0.4 | 
| R20A      |         1.3 |            -0.3 | 
| R17A      |         1.8 |            -0.2 | 
| Y21F      |         0.3 |            -0.1 | 

&lowast; G36G is a noise mutation. This is unfortunate.
The controls are as follows:

| mutations | complex_ddG | interface_ddG |
|:----------|----------:|--------------:|
| A16A |       3.4 |          -2.1 |
| G36G |      -0.5 |            -2 |
| G37G |       0   |            -0 |
| wt        |      -0   |        0 |
| G12G |       0.5 |             0 |
| K15K |      -1   |             0 |
| R17R |       0.6 |             0 |
| K46K |       0.2 |             0 | 
| T11T |       0.6 |           0.1 |       
| F45F |       0.1 |           0.1 |   
| P13P |       0.1 |           0.1 | 
| R20R |       0.1 |           0.1 |   
| Y21Y |       0.3 |           0.1 | 
| I18I |       0.9 |           0.1 | 
| I19I |       1   |           0.1 |   
| T32T |       0.8 |           0.3 |  
| V34V |       2.6 |           0.4 | 

This is likely due to the constrained minimisation.
Bar for G36G and A16A the noise is low.

## FastDesign
 
> Code: [04-fastdesign.py](04-fastdesign.py)

| mutations                                              | complex_ddG |   interface_ddG |
|:-------------------------------------------------------|------------:|----------------:|
| T11E P13W R17L Y21F T32A V34I G36A                     |         0.2 |            -6.6 |
| T11L R17L I18L I19Y R20A T32V V34I G36A                |        -2.9 |            -5.6 |
| T11L R17L I18L I19Y R20A Y21W T32V V34I G36A           |        -1.9 |            -5.4 |
| T11E P13W R17L R20V Y21W T32A V34I G36A                |         0.7 |            -5   |
| T11V R17L I18L I19F R20A Y21W T32V V34I G36A S47T      |        -2.5 |            -4.8 |
| T11S R17L I19F R20V T32V V34I G36A S47K                |          -4 |            -4.8 |
| T11L R17L I18L I19Y R20A Y21W T32V V34I G36A S47T      |        -2.9 |            -4.6 |
| T11E P13N R17L I19W R20A Y21A T32Y V34I G36A S47K      |        -4.8 |            -4.4 |
| T11E P13W R17L Y21W T32A V34I G36A S47T                |         1.5 |            -4.3 |
| T11V P13N R17L I18L I19Y R20V Y21W T32V V34I G36A S47T |          -2 |            -4.3 |
| P13N R17L R20A Y21W T32A V34I                          |        -2.6 |            -4.2 |
| T11L R17F I18L I19L R20A T32V V34I G36A S47K           |        -6.3 |            -4   |
| P13V R17L R20A Y21W T32A V34I                          |        -2.8 |            -3.7 |
| T11V R17L I18L I19F R20V Y21W T32V V34I G36A S47T      |        -0.7 |            -3.7 |
| T11V R17L I18L I19Y R20V Y21W T32V V34I G36A S47T      |          -2 |            -3.7 |
| T11E R17L I19Y R20A Y21W T32V V34I G36A S47T           |        -3.3 |            -3.6 |
| T11E R17L I19W R20A Y21A T32Y V34I G36A S47K           |        -2.4 |            -3.5 |
| T11V R17L I18L I19F Y21W T32V V34I G36A S47T           |          -3 |            -3.4 |
| R17L R20V Y21F T32A V34I                               |        -1.8 |            -3.3 |
| R17L I19W R20A Y21A T32Y V34I                          |        -3.9 |            -3.2 |
| T11L R17L I19F R20V Y21W T32V V34I G36A S47T           |        -3.6 |            -3.2 |
| T11L R17Y I18L I19L Y21W T32V V34I G36A S47T           |        -3.6 |            -3.1 |
| T11E P13N R17F I19L R20V Y21W T32V G36A S47T           |        -3.1 |            -3.1 |
| P13N R17L R20A Y21F T32A V34I                          |        -2.2 |            -3.1 |
| T11E R17L I19W R20V Y21A T32Y V34I G36A S47K           |        -4.7 |            -3   |
| T11E R17F I19L R20V T32V G36A S47T                     |        -1.8 |            -2.9 |
| T11E R17Y I19Y R20A Y21A T32Y V34I G36A S47H           |        -3.1 |            -2.9 |
| Y21A Y23W G36A R39D A48L T54V G56A                     |       -10.7 |            -2.9 |
| T11V R17L I19W R20V Y21A T32Y V34I G36A S47K           |        -3.7 |            -2.8 |
| R17Y R20V Y21F T32A                                    |          -1 |            -2.8 |
| T11V R17F I18L I19L R20A Y21W T32V G36A                |          -2 |            -2.8 |
| I18L                                                   |         1.4 |            -2.6 |
| R17L I19W R20V Y21A T32Y V34I S47K                     |        -4.9 |            -2.5 |
| R17L I19W R20V Y21W T32V V34I S47T                     |        -1.8 |            -2.5 |
| Y21A Y23F G36A R39D A48L T54V G56S                     |       -11.9 |            -2.5 |
| T11S R17L Y21W T32A V34I G36A                          |        -1.7 |            -2.4 |
| T11V R17L I19L R20V Y21W T32V V34I G36A S47T           |        -2.5 |            -2.4 |
| T11V R17L I19L R20V Y21W T32V V34I G36A                |        -1.9 |            -2.4 |
| T11E R17F R20V Y21W T32A V34I G36A                     |        -0.1 |            -2.3 |
| P13V R17L R20A Y21F T32A V34I                          |        -2.2 |            -2.2 |
| P13V R17Y R20A Y21W T32A                               |        -1.5 |            -2.2 |
| R17Y I18L Y21F T32A                                    |        -1.3 |            -2.2 |
| T11E P13N R17F I18L R20Y Y21W T32P V34I G36A           |        -5.1 |            -2.2 |
| T11L R17Y I18L I19V R20A G36A                          |        -2.7 |            -2.2 |
| T11E R17F R20A Y21W T32A V34I G36A                     |        -1.7 |            -2.1 |
| T11L R17L I18L Y21W T32A V34I G36A                     |        -0.4 |            -2.1 |
| T11E R17L R20A Y21F T32A V34I G36A                     |        -0.7 |            -2.1 |
| T11E K15R R17L I19W R20V Y21A T32Y V34I G36A S47K      |        -2.7 |            -2.1 |
| I19W                                                   |         1.6 |            -2.1 |
| P13V R17F R20A Y21W T32A V34I                          |        -1.7 |            -2.1 |
| T11E R17L R20V Y21F T32A V34I G36A                     |          -0 |            -2   |
|                                                        |        -0.5 |            -2   |
| T11E R17F R20V Y21W T32A V34I G36A S47T                |           1 |            -2   |
| R17L R20A Y21F T32A V34I                               |          -2 |            -2   |
| T11E R17L I18L R20V Y21A T32Y V34I G36A S47T           |        -0.7 |            -2   |
| R17Y Y21F T32A                                         |        -2.2 |            -2   |
| T11E P13N R17F R20V Y21W T32A V34I G36A S47T           |         0.6 |            -1.9 |
| P13V R17L R20V Y21F T32A V34I                          |        -1.7 |            -1.9 |
| P13N R17L I18L Y21W T32A V34I                          |        -1.2 |            -1.9 |
| R17L I18L Y21F T32A V34I                               |        -0.8 |            -1.8 |
| T11V R17Y I18L I19L Y21W T32V V34I G36A                |        -3.3 |            -1.8 |
| R17F R20V Y21W T32A V34I S47T                          |        -1.6 |            -1.8 |
| T11L R17L R20A Y21W T32A V34I G36A S47T                |          -1 |            -1.8 |
| T11V R17F I19L R20V Y21W T32V G36A S47T                |          -3 |            -1.8 |
| T11V R17F R20V Y21W T32A V34I G36A                     |           1 |            -1.8 |
| T11V R17F I18L I19L Y21W T32V V34I G36A S47T           |        -3.2 |            -1.8 |
| T11V R17Y I18L I19L R20V Y21W T32V G36A S47T           |        -1.3 |            -1.8 |
| P13N R17L R20A Y21A T32M V34I                          |         1.4 |            -1.8 |
| T11V R17F I18L I19L Y21W T32V G36A S47T                |        -2.6 |            -1.7 |
| T11S R17F R20V Y21W T32A V34I G36A                     |        -0.6 |            -1.7 |
| P13N R17L R20V Y21W T32A V34I                          |        -1.9 |            -1.7 |
| T11L R17L R20A T32A V34I G36A                          |        -1.5 |            -1.7 |
| T11S R17F R20A Y21W T32A V34I G36A                     |          -2 |            -1.7 |
| P13D R17Y R20A Y21W T32A                               |        -2.3 |            -1.7 |
| T11E P13D R17L R20A Y21F T32A V34I G36A                |        -0.7 |            -1.7 |
| R17L I18L R20A T32A V34I                               |        -1.6 |            -1.7 |
| T11E R17Y Y21F T32A G36A                               |        -0.5 |            -1.7 |
| R17L R20A T32A V34I                                    |        -0.9 |            -1.7 |
| R17F I18L Y21W T32A V34I                               |        -2.1 |            -1.6 |
| R17F I18L I19L Y21W T32V V34I G36A                     |        -2.9 |            -1.6 |
| R17F I18L Y21W T32A V34I S47T                          |        -1.6 |            -1.6 |
| P13D R17L I18L Y21F T32A V34I                          |          -2 |            -1.6 |
| T11S R17F Y21W T32A V34I G36A                          |        -1.2 |            -1.6 |
| T11V R17L R20V Y21W T32A V34I G36A S47T                |        -0.5 |            -1.6 |
| R17L I18L Y21F T32A V34I S47T                          |        -1.3 |            -1.5 |
| T11S R17L I19W R20V Y21A T32Y V34I G36A S47T           |        -4.3 |            -1.5 |
| R17L R20A Y21W T32A V34I                               |        -2.6 |            -1.5 |
| T11E R17L I18L Y21F T32A V34I G36A                     |         0.2 |            -1.5 |
| P13E                                                   |         1.7 |            -1.5 |
| R17L I18L T32A V34I S47T                               |        -1.2 |            -1.5 |
| T11E P13V R17F I19L R20V Y21W T32V G36S S47T           |        -1.7 |            -1.5 |
| T11E R17F I18L T32A V34I G36A                          |        -1.3 |            -1.5 |
| R17L I18L T32A V34I                                    |        -1.5 |            -1.4 |
| T11V R17L I18L Y21W T32A V34I G36A S47T                |          -1 |            -1.4 |
| R17L Y21F T32A V34I                                    |        -2.4 |            -1.4 |
| R17L I18L Y21W T32A V34I G36A                          |        -1.4 |            -1.4 |
| R17M                                                   |         1.1 |            -1.4 |
| R17L I19W R20V Y21W T32A V34I S47T                     |        -2.3 |            -1.4 |
| T11E R17Y I18L Y21W T32A V34I G36A S47T                |        -1.1 |            -1.3 |
| T11V A16V R17L I18L R20V Y21A T32Y V34I G36A S47K      |        -1.8 |            -1.3 |

For analysis of repetitions and epistasis, see email.
 
## Future

Possible improvements:

* I did not run an MD run or Rosetta Backrub to get alternative conformations of the protein for better sampling.
* I forgot to calculate the monomer ∆G...
* I can re-run everything on the variant without cysteines...