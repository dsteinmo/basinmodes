# basinmodes

This MATLAB/GNU Octave code repository can be used to calculate basin-scale free oscillations in natural bodies of water, in particular their frequencies and spatial structure as the solution ("modes") of an eigenvalue problem. Visualization of the oscillations is also facilitated by these scripts. See below for details.

## Getting Started

    $ git clone https://github.com/dsteinmo/basinmodes.git
    $ cd basinmodes
    $ matlab -nodesktop -nosplash          # OR: octave-cli

From here (In MATLAB/Octave), cd into the various test case folders and run the driver scripts. The driver scripts have the same name as the folders.

    > cd csanady_internalmodes
    > csanady_internalmodes

All physical and computational parameters (and mesh) are specified within the driver scripts. Implementation details are in basinmodes_curved.m (or basinmodes.m) that lives in the root folder.

The circular-basin test cases (Csanady; Lamb) use the circlemesh_nohole.mat mesh that is checked into the repository root folder. For Lake Ontario test case, you will need to generate the mesh file and copy it into the ontario_surfacemodes/ folder.

For machines with gmsh installed, this can be done (in bash) vioa:

    $ chmod +x geo2msh.sh
    $ ./geo2msh.sh ontario.geo 

Then (in Matlab/Octave), run `msh2mat`. Note: You will need to update the top of `msh2mat.m` to match your input and output filenames.

    > msh2mat

Then copy the resulting .mat file into the ontario_surfacemodes/ folder and ensure the driver script (e.g., ontario_surfaces.m) is referencing the .mat file appropriately.

## Visualization

To interactively scroll through the modes (with `⇐` and `⇒` arrow keys) run the `scrollmodes.m` script:

    > scrollmodes

* `animate_mode.m` can be modified and used to animate one period of a particular mode that is specified at the top of the script.

## Citing this work

Results from using the code herein was published in 'Ocean Nodelling' in 2019. BibTeX:
```
@article{steinmoeller2019calculating,
  title={Calculating basin-scale free oscillations in lakes on a rotating Earth},
  author={Steinmoeller, DT and Stastna, M and Lamb, KG},
  journal={Ocean Modelling},
  volume={139},
  pages={101403},
  year={2019},
  publisher={Elsevier}
}
```