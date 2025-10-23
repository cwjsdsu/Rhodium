Code RHODIUM is a post-processor for BIGSTICK. Eventually it will compute
generate densities (rhos) hence the name. Right now all it does is 
dot products. It allows for different initial (_i) and final (_f) bases.
The bases can be the same, but need to be read in nonetheless.

Because I am modifying the BIGSTICK routines, there are a number of BIGSTICK files here, 
most of them unused so far. These are all of the form bXXXXXX.f90, e.g., bbasislib1.f90,
breorthog.f90, etc.. You can ignore them; they are here for further development by me.

To compile, 

make gfortran

which will create the executable 

rhodium.x

As of this version, only serial is working.

To run RHODIUM, you must first use BIGSTICK to create both .wfn files and a .bas file.

To understand the basis, the key file is 

rhbasislib.f90

Currently the basis is stored in bit representation, in multiple words. 

The proton SDs are stored in psdlist_i and neutron in nsdlist_i (for the initial basis)
These are stored as bits. To convert to a list of occupied orbits, you can use the routine

convert_bitrepsd2occ


There are nxsd_i(it) proton/neutron slater determinants, it =1 (proton), 2 (neutron)

To map from the proton and neutron slater determinant indices to the final basis used by Bigstick, 
let ip be the proton label and jn the neutron label. Then

basis index = pstart_i(ip) + nstart_i(jn)

TO RUN you need both a .wfn (from a normal run in Bigstick) and a .bas file (option 'b') from the main menu. 
Actually, you need these for both the initial and final states. 



