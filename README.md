#Getting Started

    git clone https://github.com/dsteinmo/basinmodes.git
    cd basinmodes
    matlab -nodesktop -nosplash

From here (In Matlab/Octave), cd into the various test case folders and run the driver scripts. The driver scripts have the same name as the folders.

    cd csanady_internalmodes
    csanady_internalmodes

All physical and computational parameters (and mesh) are specified within the driver scripts. Implementation details are in basinmodes_curved.m (or basinmodes.m) that lives in the root folder.

The circular-basin test cases (Csanady & Lamb) use the circlemesh_nohole.mat mesh that is checked into the repo. For Lake Ontario test case, you will need to generate the mesh file and copy it into the ontario_surfacemodes/ folder.

For machines with gmsh installed, this can be done with (in bash):

    . geo2msh.sh ontario.geo 

Then (in Matlab/Octave):

    msh2mat

Then copy the resulting .mat file into the ontario_surfacemodes/ folder and ensure the driver script is referencing the .mat file appropriately.

#Visualization

To interactively scroll through the modes (with left and right arrow keys) run the scrollmodes.m script:
	
    scrollmodes

animate_mode.m can be modified and used to animate one period of a given mode.

