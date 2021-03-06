# Step A: extraction of the tetrahedral representation of the channel


### Before you start
The identification of the (tetrahedral) connected components is integrated within [NanoShaper](https://electrostaticszone.eu/downloads). While the packages and the procedure required to run this version are the virtually same as NanoShaper 0.7, this version can be installed in more recent operating systems (tested on MacOS Big Sur and Monterey and on Ubuntu 20.04) and with more recent packages (e.g., CGAL 5.3.1). A short guide of how to set up NanoShaper in Ubuntu 20.04 is included. A similar procedure can be applied, for example, to install NanoShaper in macOS.

### Example files
Two conformations can be found in `example/IonicChannels/mscl/5ligand/`:
- The file `NanoShaper_input.xyzr` contains the atomic centers and radii, one per line.
- The file `conf.prm` is the configuration file that provides NanoShaper with its input (e.g., type of molecular surface to adopt and probe radius).
These two files are what NanoShaper expects to run properly.

### On compilation and run
Note that, to make compilation easier, a CMakeLists file and a shell file (`Chanalyzer_cmake.sh`) are included in the folder. An example of shell file to run the exe file generated by compilation is `Chanalyzer_example.sh`. The following flags can be set:
- `-conf`: the name of the configuration file.
- `-find_geo_pockets`: if you want to find candidate geometric pockets (when set to 1) or just run NanoShaper (when set to 0).
- `-num_c2r`: number of candidate pockets to save, from the largest (in terms of volume) to the smallest ones.
 
### Output files
 In addition to the usual output generated by NanoShaper, the following files are expected to be generated:
 - `ED/CCs/*.OFF`: N largest (w.r.t. volume) tetrahedral connected components
 - `ED/ED_tet/*.txt`: for each tetrahedral connected component in `ED/CCs`, there is a corresponding text file that contains, for each line, four triplets: they are the four points that generate tetrahedra of such connected component. 
