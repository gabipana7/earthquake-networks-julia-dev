## earthquake-networks-julia-dev

## DEVELOPMENT 
### This project aims to describe Seismic Regions using elements of Complex Networks
### This is an enhanced Julia port of the original project done in Python

---
# Exploratory Data Analysis
For an extensive review of the data collection and cleaning process, as well as graphical exploration of the data, codes and figures are presented in the complementary project: https://github.com/gabipana7/seismic-exploratory-data-analysis

---
# Network Generation
First of all networks are generated and their connectivity distribution is analyzed.

This is done with *networks_parameter_dependency.jl*

---
# Netwoks and Motifs Generation
With *networks_generator.jl* networks are generated for the best cube sizes, obtained from the parameter dependency.

With *motifs_generator.jl* triangle and tetrahedron motifs are generated for these networks, with a cutoff for micro-earthquakes (data is trimmed to satisfy that magnitude > 2). 

This computation is done with a slightly adapted Python version of the Nemomap software (https://github.com/zicanl/NemoMapPy)


---
# Motifs Analysis
With *motifs_analysis.jl* the motifs are analyzed as such:

Triangles
- the area of each individual motif is calculated
- for each individual motif, the total (and mean) energy released by the earthquakes that occured in each of its nodes is computed
- for each individual motif, the area is weighted by its total(mean) energy
- the distribution of the results is computed and analyzed with the Python Powerlaw package

Tetrahedrons
- the volume of each individual motif is calculated
- for each individual motif, the total (and mean) energy released by the earthquakes that occured in each of its nodes is computed
- for each individual motif, the volume is weighted by its total(mean) energy
- the distribution of the results is computed and analyzed with the Python Powerlaw package
