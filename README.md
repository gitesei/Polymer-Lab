[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/mlund/particletracking/master)

# Polymer-Lab

In this computer lab, you will invetigate polymer mediated interactions between non-adsorbing (hard and inert) surfaces. Two main scenarios will be studied:

<ol type="a">
  <li><b>Net interactions in a solution containing free (ungrafted) polymers.<\b>
We will here limit ourselves to the case of theta solvents, and utilize the standard ideal polymer model, in which the monomers are point-like, i.e they do not interact with each other (but they are still excluded from the surfaces). Run calculations with various choices for the polymer length. Plot monomer density distributions at some chosen separations, and comment upon the results. Also compare the free energy, and net pressure curves that you obtain for different polymer lengths. Limit the degree of polymerization ($r$) to $3 < r < 1000$. Also limit the investigated separations ($h$) to $3 < h < 500$ (where we use the bond length as length unit). Also plot a few free energy curves, with various polymer lengths, but use radius of gyration, rather than bond length, on the $x$-axis. Discuss the results. Note that the the bulk monomer concentration is fixed (by the code). How would the free energy curves change if you had (say) doubled the bulk density? 
0.  Specifically plot a free energy curve for 601-mers, and verify that you obtain the same results (with bond lengths on the $x$-axis) as in the "polydispersity comparison graph" that you have been given. Discuss and interpret these results, i.e. how and why the interactions differ when the solution contains polydisperse polymers.
0.  Using the Derjaguin Approximation, the net force, $F$, between (large) spherical particles of radius R, can be etimated from the net free energy per unit area between flat surfaces, $\Delta g_S$, as $F = R\pi\Delta g_S$. How would you obtain the net interaction free energy, $W$, between such particles? If we neglect other (non polymer-mediated) interactions, make a qualitative comparison between the second virial coefficients ($B_{22}$) for particles in monodisperse and polydisperse polymer solutions. In other words, in which case would you expect a more attractive $B_{22}$?</li>
  <li><b>Interactions between polymer grafted surfaces.i<\b>
Also in this case, the polymers are assumed to be ideal. However, you will also receive results obtained with a good solvent, as modelled by hard-sphere monomers. 
Use the same "$h$ and $r$" limitations as above ($3 < r < 1000$; $3 < h < 500$).<br>
Plot density profiles, free energy curves, and pressure curves, for some chosen polymer lengths and grafting densities. Discuss the results.<br>
Specifically plot a free energy curve for 100-mers, with a reduced grafting density of 0.1, and verify that you obtain the same results as in the "solvent comparison interaction free energy graph" that you have been given.<br>
Also plot the density profile, at h = 100 (bond lengths) and compare with the provided "solvent comparison density profile graph".<br>
Sketch, directly in that graph, how the density profile in a good solvent would look (roughly) had the surface grafting density been half as high. How would you expect that interaction free energy curve would change if the surface were attractive to the monomers?</li>
</ol>
 
## Prerequisites

- No prior knowledge in Python is required, but familiarity with programming concepts is helpful.
- A laptop connected to the internet (eduroam, for example) and running Unix, MacOS, or Windows and with Anaconda installed, see below.

If you have little experience with Python or shell programming, the following two tutorials may be helpful:

- https://swcarpentry.github.io/shell-novice
- https://swcarpentry.github.io/python-novice-inflammation

## Preparation before the lab

1. Install [miniconda3](https://conda.io/miniconda.html).
2. [Download](https://github.com/mlund/particletracking/archive/master.zip) the lab material
   (this github repository) and unzip.
3. Install and activate the `particletracking` environment described by the file [`environment.yml`](/environment.yml)
   by running the following in a terminal:

   ```bash
   conda env create -f environment.yml
   source activate particletracking
   ```
N.B. Not available for Windows. 

[Further Information](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

## How to launch the notebooks

~~~ bash
jupyter-notebook A/A.ipynb
jupyter-notebook B/B.ipynb
~~~
