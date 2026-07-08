# Multi-Phase-Field-for-FRC-using-Puck-theory

The codes presented here are used in the paper "Revisiting Multi-Phase field for FRCs using Puck failure theory". The codes can be used with AT2, AT1 models. We have tested the files for version Abaqus2020, Abaqus2022, Abaqus2024. 
All the material properties required are presented in the pfall.f file, and the input file corresponding to each example in the paper are given. 

Please note that if your initial mesh is Bigger than 70000, Third line in the module Kvisual "UserVar(70000,16,4)" has to be changed to the upper limit of the number of the mesh. 

If you are using this code for the academic research or industrial purpose, please cite our paper 



1. Pavan Kumar Asur, Rafeal Fleischhaker, Aamir Dean, Heinz E Pettermann "Revisiting Multi-phase field model for FRCs using Puck failure theory", Under Revision to Composites Structures.
2. A. Dean, Pavan Kumar AV, J. Reinoso, C. Gerendt, E.Mahdi, M. Paggi, R. Rolfes, "A multi-phase field fracture model for long fibre reinforced composites based on the Puck theory of failure" composite structures Volume 251, 1 November 2020, 112446.
3.  Pavan Kumar AV, A. Dean, J. Reinoso, M. Paggi, "A Multi Phase-Field-Cohesive Zone Model for Laminated Composites: Application to Delamination Migration" composite structures,Volume 276, 15 November 2021, 114471.


Authors: 

Pavan Kumar Asur: "Pavan-asur", pavan.kumar@tuwien.ac.at, asurpavankumar@gmail.com

Aamir Dean: "aamir-dean", a.dean@sustech.edu, dr.aamir.dean@gmail.com
