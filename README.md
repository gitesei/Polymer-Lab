[![Binder](https://mybinder.org/badge.svg)](https://mybinder.org/v2/gh/gitesei/Polymer-Lab/master)

# Polymer-Lab

In this computer lab, you will invetigate polymer mediated interactions between non-adsorbing (hard and inert) surfaces. Two main scenarios will be studied:

## A. Net interactions in a solution containing free (ungrafted) polymers.
We will here limit ourselves to the case of theta solvents, and utilize the standard ideal polymer model, in which the monomers are point-like, i.e they do not interact with each other (but they are still excluded from the surfaces).<br> 
Further instructions in `A/A.ipynb`.

## B. Interactions between polymer grafted surfaces.
Also in this case, the polymers are assumed to be ideal. However, you will also receive results obtained with a good solvent, as modelled by hard-sphere monomers.<br> 
Further instructions in `B/B.ipynb`.

Author: Jan Forsman<br>
E-mail: jan.forsman@teokem.lu.se
 
### Prerequisites

- No prior knowledge in Python is required, but familiarity with programming concepts is helpful.
- A laptop connected to the internet (eduroam, for example) and running Unix, MacOS, or Windows and with Anaconda installed, see below.

If you have little experience with Python or shell programming, the following two tutorials may be helpful:

- https://swcarpentry.github.io/shell-novice
- https://swcarpentry.github.io/python-novice-inflammation

### Preparation before the lab

1. Install [miniconda3](https://conda.io/miniconda.html).
2. [Download](https://github.com/mlund/particletracking/archive/master.zip) the lab material
   (this github repository) and unzip.
3. Install and activate the `polymerlab` environment described by the file [`environment.yml`](/environment.yml)
   by running the following in a terminal:

   ```bash
   conda env create -f environment.yml
   source activate polymerlab
   ```
N.B. Not available for Windows. 

[Further Information](https://conda.io/docs/user-guide/tasks/manage-environments.html#creating-an-environment-from-an-environment-yml-file)

### How to launch the notebooks

~~~ bash
jupyter-notebook A/A.ipynb
jupyter-notebook B/B.ipynb
~~~
