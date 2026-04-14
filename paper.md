---
title: 'DynDom1D_Python: An Open-Source Python tool for the analysis of domain movements in proteins'
tags:
  - Python
  - protein
  - dynamics
  - hinge-bending
  - hinge-axis
authors:
  - name: Jien Lim
    affiliation: 1
  - name: Hugh Millard
    orcid: 0009-0004-8963-1737
    affiliation: 1
  - name: Steven Hayward
    orcid: 0000-0001-6959-2604
    corresponding: true
    affiliation: 1
affiliations:
 - name: School of Computing Sciences, University of East Anglia, Norwich NR4 7TJ, U.K.
   index: 1

date: 1 August 2025
bibliography: paper.bib
---

# Summary

Here we introduce, DynDom1D_Python, an open-source Python implementation based on the original DynDom Fortran program for the analysis of domain movements in proteins. DynDom works on a single protein chain and can be used when two structures of the same protein are available representing a conformational change. If appropriate it describes the conformational change in terms of the relative rotation of quasi-rigid regions called "dynamic domains" by way of hinge axes (more precisely, interdomain screw axes) and hinge-bending residues (more precisely, interdomain bending residues). This new implementation improves on the original standalone version in a number of ways.

# Statement of need

Protein function normally involves interaction with other molecules, this interaction often effecting a change in conformation. Proteins are macromolecules that usually comprise many tens of thousands of atoms and conformational change occurs over a vast spatial scale. On the global scale is the domain movement where each domain comprises a significant portion of the whole protein. Protein domain movements form an important class of movements in proteins, perhaps the most important, and they are often engaged in all aspects of protein function, e.g. an enzyme trapping of a small substrate in its interdomain cleft. In principle, understanding how domain movements are controlled during function may be useful for drug design, e.g. by designing a drug to bind to a hinge-bending region in order to disrupt a functional hinge-bending movement.

The original Fortran implementation of the DynDom program [@Hayward:1998;@Lee:2002] is available to download and run at the DynDom webserver [@Lee:2003] and to download from Collaborative Computational Project 4 (CCP4) [@Agirre:2023], where it has been used by protein X-ray crystallographers since 1999. Although versions of DynDom for the analysis of multimeric proteins (comprising multiple chains) have been released, namely DynDom3D [@Poornam:2009] and DynDom6D [@Veevers:2019], it is the 1D version (1D indicating that it analyses a single protein chain) that remains the most popular. An updated, open-source version of DynDom1D with a less restrictive licence and implemented with a widely used programming language, will enable others to contribute.

# State of the field

There are two main approaches to the task of determining dynamic domains and hinge axes from a pair of structures. The DynDom tools and the DomainFinder tool [@Hinsen:1999] cluster parameters relating to rigid body motions to determine dynamic domains. The HingeFind approach of Wriggers and Schulten [@Wriggers:1997] determines domains using an adaptive best-fitting procedure and is run from within X-PLOR or the VMD environment. There are distinct differences between all three tools. DomainFinder [@DomainFinder] performs normal mode analysis on a single structure to determine its modes of movement and does not directly accept two structure files. HingeFind determines the rotation axis that best describes the domain movement, whereas DynDom and DomainFinder determine the screw axis. DynDom1D is unique in being able to determine hinge-bending residues as part of an overall coherent process. As stated above hinge-bending regions are potential binding sites for drug molecules.


# Research Impact Statement

The Fortran version of DynDom1D was released in 1998 and as part of CCP4 from 1999. The first DynDom article [@Hayward:1998] has been cited by more than 950 articles (Google Scholar) and a subsequent article that improved and extended some of the basic methodologies, more than 300 times. The number of citations represents a minimum, as being part of CCP4, DynDom has been used by X-ray crystallographers with citation to CCP4 only. The majority of citations are from research articles by structural biologists who have solved the structure of their target protein in a functional state that exhibits a difference in conformation to one already known. DynDom is used to characterise this difference in terms of domains moving via hinge-bending. Target systems are often those critically involved in health conditions, e.g. neurological conditions [@Gangwar2025] and cancer [@Sotomayor2012;@HDAC62025].

# Method

The basic idea behind the DynDom approach is that domains can be recognised as regions that move as quasi-rigid bodies. DynDom1D reads two protein chain structures that represent a conformational change. These structures can be experimentally determined using methods such as X-ray crystallography, Nuclear Magnetic Resonance Spectroscopy, or cryo-electron microscopy. The two structures can also represent a set of atomic displacements derived from computational methods such as normal mode analysis or molecular dynamics simulation[@Hayward_Retro2023]. DynDom is a coarse-graining method that aims to model, when appropriate, atomic displacements in terms of the relative screw movements of a small number of "dynamic domains" moving as quasi-rigid regions, their relative movement being facilitated and controlled by flexible hinge bending regions. The process comprises the following steps:

1. Read in the two structures as specified by the PDB accession codes and chain identifiers.
2. Superpose the two structures using least-squares best fitting on the main-chain atoms N, CA, and C.
3. Use a sliding window to generate overlapping main-chain segments (default 5 residues).
4. Determine the "rotation vector" of each main-chain segment, plotting components as coordinates in the "rotation space".
5. Use K-means clustering to determine clusters of rotation points indicating regions possibly rotating within the same rigid body.
6. Use connected set algorithm to determine if regions from each cluster are spatially connected.
7. Split regions from each cluster into dynamic domains (must have at least as many residues as the minimum domain size - default 20 residues) if they are not spatially connected.
8. Check that that ratio of external to internal displacement of all dynamic domains that are connected through the main-chain is greater than the "ratio" threshold (default, 1.0).
9. If the ratio test is passed for all, determine screw axes.
10. Determine hinge bending regions.
11. Assign mechanical hinges and percentage "closure".
12. Output results including PDB files for display with PyMol using a PyMol script file (.pml file).

Further details on the basic methodology can be found in the original DynDom publications[@Hayward:1998;@Lee:2002]. An important feature is its ability to assign "mechanical hinges". These are hinge-bending regions that are close to (within 5.5 A) the hinge axis [@Hayward:1999]. In order to properly control a domain movement, two or more mechanical hinges are separated in space determining the direction and location of the screw axis. Often such arrangements create a "closure" motion whereby the interdomain cleft narrows or widens depending on the direction of rotation. DynDom1D_Python assigns mechanical hinges, and the percentage closure motion (0% closure implies 100% twisting) and also gives a point on the screw axis and its direction in terms of the unit vector.  

# Software design

 The aim of this new build was to create a Python implementation of the Fortran version of DynDom1D, to simplify and update the input and output, to incorporate methodological improvements from DynDom3D and DynDom6D, and to provide greater interoperability by using tools from popular packages. In terms of the user experience, it has a much simplified input (requiring only the protein data bank (PDB) accession codes and the chain identifiers of the two structures), and an output that produces a script file for display and animation of the domain movement using PyMol[@pymol], the most widely-used molecular graphics tool - the Fortran version produces scripts for the now outdated RasMol software. Tools for reading PDB files and structure superposition are from GEMMI [@Wojdyr:2022], and kmeans clustering using the scikit-learn tool[@pedregosa2011scikit:2011]. To facilitate readability and maintainance, the code has classes containing associated attributes and methods. This new implementation uses the language of choice of the biological sciences community and popular packages, making it maintainable and future-proof. 

# Example Result - Adenylate Kinase

\autoref{fig:Fig1} shows the results of the domain movements that occur in the protein adenylate kinase upon binding of the small molecules ADP and AMP. It shows how the protein closes up upon these two molecules. Details regarding these movements are given in the "info" file where it states that the rotations are 53 degrees for the red moving domain, and 46 degrees for the yellow moving domain relative to the central blue fixed domain. Both motions are controlled by a pair of separated mechanical hinges and are closure motions (99% and 97%, respectively).


![DynDom1D_Python result visualised using PyMol, for the movement between 4AKE(A) and 2ECK(B), the free and ADP,AMP bound structures of adenylate kinase. The "fixed" domain in the centre is coloured blue and the two "moving" domains are coloured red and yellow. The hinge bending regions are coloured green (DynDom always colours bending regions green) and the interdomain screw axes are depicted as arrows, with the shaft coloured the same as the fixed domain and the arrow head colour indicating the associated moving domain.\label{fig:Fig1}](adenylate_kinase.png)

# AI usage disclosure

All methodologies were conceived and designed by the authors. Claude Sonnet 4.0 was used to assist in translating some legacy Fortran functions into Python. All AI generated code was reviewed, checked and edited by the authors. Testing of the output was performed by the authors. AI was not used in the preparation of this article.

# References

