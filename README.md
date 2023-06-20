# Automatic Johnson Cluster General Algorithm

A simple algorithm to create chemical structures with the Jhonson Cluster conformations.
Please see [Article](none.com).

# Installation

Simply download .py file from the Github: [Link](https://github.com/Rufox/Automatic-Johnson-Cluster-General-Algorithm)

# Requirements

 - Python 3.8 (tested)
	 - Sympy library
	 - Numpy library
	 - Configparser library
-	Files
	-	Parameters file
# Usage
To initialise the program use the following command line

    python AJCG.py Config.in

The algorithm will create the following files:

 - **Results_Combinations.xyz** : **XYZ** file with all the possible conformations created using the **Config** file parameters. Each structure identified with a special tag.
 - **InputsGaussian/**
	 - *PermutXX.com* : Gaussian input file for an specific structure. **XX** references to the special laber in the **XYZ** file

# Configuration File
The configuration file contains all the necessary parameters for the construction of a Jhonson clsuter, as well as the chemical equation with which it will be populated.

## GENERAL

 - **Shape**:  parameter to specify the shape of a specific Jhonson cluster. Clusters are constructed using rings called "floors". They are then arranged vertically to form the final structure. Each ring is represented by the number of vertices (a square would be **4**). To rotate a ring a negative label is given (**-4**). Each floor is separated by a comma (**,**). For example: 
	 - [Icosahedron](https://en.wikipedia.org/wiki/Icosahedron "Icosahedron"): **1, 5,-5,1**
	 - [Triangular bipyramid](https://en.wikipedia.org/wiki/Triangular_bipyramid "Triangular bipyramid"): **1, 3, 1**.
	 - [Gyrobifastigium](https://en.wikipedia.org/wiki/Gyrobifastigium "Gyrobifastigium"): **2, 4, -2**
 - **Distances**: parameter to control the bond distance between the vertices for each floor (in Angstrom). This parameter is directly correlating with "shape" and MUST be the same length. For example, to build a pyramid of 4 floors, the parameter would be:
	 - shape = **3, 3, 3, 1**
	 - distances = **4, 3, 2, 1** 
 - **Height**: Parameter to control the bond distances between each floor. Of the length : length(shape) -1. Very high values will elongate the structure.

## Cluster

 - **Chemical_formula**: parameter to specify the chemical structure to poblate the jhonson cluster prevously created. The format is **X N**, where **X** is any atomic symbol and **N** the number of that atom. Water (H2O) would be **H 2 O 1**

## Software

 - **software**: Chemical software for an input to be created. Gaussian, ORCA.
 - **core**: Number of CPUs used in the calculation for the software input
 - **memory**: RAM usage (GB) for software input.
 - **charge_multi**: Charge and Multiplicity of the molecule.
 - **header**: Specific order to be given to the software (optimization, scf calculation, etc)