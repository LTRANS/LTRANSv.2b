# LTRANSv.2b
**Larval TRANSport Lagrangian model is an off-line particle-tracking model that is forced with Regional Ocean Modeling System (ROMS) results.**

The code was originally released here:
http://northweb.hpl.umces.edu/LTRANS.htm

Code in the LTRANS GitHub organization is based off of code developed by Elizabeth North, Zachary Schlag, and Ian Mitchell with input from Chris Sherwood and Scott Peckham.  

-------------------------------
#### The following papers and User's Guide should be cited:  

North, E. W., E. E. Adams, Z. Schlag, C. R. Sherwood, R. He, S. Socolofsky. 2011. Simulating oil droplet dispersal from the Deepwater Horizon spill with a Lagrangian approach. AGU Book Series: Monitoring and Modeling the Deepwater Horizon Oil Spill: A Record Breaking Enterprise  [link to paper]( http://onlinelibrary.wiley.com/doi/10.1029/2011GM001102/summary)

North, E. W., Z. Schlag, R. R. Hood, M. Li, L. Zhong, T. Gross, and V. S. Kennedy. 2008. Vertical swimming behavior influences the dispersal of simulated oyster larvae in a coupled particle-tracking and hydrodynamic model of Chesapeake Bay. Marine Ecology Progress Series 359: 99-115  [link to paper](https://doi.org/10.3354/meps07317)

Schlag, Z. R., and E. W. North. 2012. Lagrangian TRANSport model (LTRANS v.2) User’s Guide. University of Maryland Center for Environmental Science, Horn Point Laboratory. Cambridge, MD. 183 pp.  [link to guide](https://github.com/LTRANS/LTRANSv.2b/blob/master/LTRANSv2_UsersGuide_6Jan12.pdf)

-------------------------------
#### External dependencies and programs:

LTRANS v.2b requires NetCDF libraries and uses the following programs to calculate random numbers (Mersenne Twister) and fit tension splines (TSPACK). Because LTRANS v.2 reads-in ROMS-generated NetCDF (.nc) files, it requires that the appropriate NetCDF libraries be installed on your computer (see files and links below). Also, please note that although the Mersenne Twister and TSPACK programs are included in the LTRANS v.2b in the Random_module.f90 and Tension_module.f90, respectively, they do not share the same license file as LTRANS v.2b. Please review and respect their permissions (links and instructions provided below).

**Windows Visual Fortran NetCDF libraries:**
These NetCDF files that are compatible with Visual Fortran were downloaded from the Unidata NetCDF Binaries Website for LTRANS v.1. The NetCDF 90 files were downloaded from Building the F90 API for Windows for the Intel ifort compiler website. The VF-NetCDF.zip folder contains README.txt that describes where to place the enclosed files. If these files do not work, you may have to download updated versions or build your own by following the instructions at the UCAR Unidata NetCDF website.

**Linux NetCDF libraries:**
Linux NetCDF libraries. Linux users will likely have to build their own Fortran 90 libraries using the source code/binaries that are available on the UCAR Unidata NetCDF website.

**Mersenne Twister random number generator:**
This program was recoded into F90 and included in the Random_module.f90 in LTRANS. See the Mersenne Twister Home Page for more information about this open source program. If you plan to use this program in LTRANS, please send an email to: m-mat @ math.sci.hiroshima-u.ac.jp (remove space) to inform the developers as a courtesy.

**TSPACK: tension spline curve-fitting package:**
This program (ACM TOMS Algorithm 716) was created by Robert J. Renka and is used in LTRANS as part of the water column profile interpolation technique. The original TSPACK code can be found at the link to the left and is copyrighted by the Association for Computing Machinery (ACM). With the permission of Dr. Renka and ACM, TSPACK was modified for use in LTRANS by removing unused code and call variables and updating it to Fortran 90. The modified version of TSPACK is included in the LTRANS source code in the Tension Spline Module (tension_module.f90). If you would like to use LTRANS with the modified TSPACK software, please read and respect the ACM Software Copyright and License Agreement. For noncommercial use, ACM grants "a royalty-free, nonexclusive right to execute, copy, modify and distribute both the binary and source code solely for academic, research and other similar noncommercial uses" subject to the conditions noted in the license agreement. Note that if you plan commercial use of LTRANS with the modified TSPACK software, you must contact ACM at permissions@acm.org to arrange an appropriate license. It may require payment of a license fee for commerical use.

-------------------------------
#### Example input files

These files can be used to test LTRANS v.2b. They include examples of particle location and habitat polygon input files (.csv) and ROMS grid and history files (.nc) that are needed to run LTRANS v.2b. Many thanks to Wen Long for sharing the ROMS .nc files. The LTRANS v.2b code above is configured to run with these input files. Note: please download the tar (LTRANSv2.tgz) history files (clippped_macroms_his_*.nc) files between the hours of 5 pm and 6 am Eastern Standard Time because of their large size.  These files are too large to be included in the GitHub repository.  

**Initial particle locations and habitat polygon input files needed for initializing LTRANS. The habitat polygon input files (End_polygons.csv and End_holes.csv) are only required when the settlement module is turned on.**

[Initial_particle_locations.csv](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/Initial_particle_locations.csv)  
[End_polygons.csv](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/End_polygons.csv)  
[End_holes.csv](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/End_holes.csv)


**ROMS grid (for LTRANS v.2b input). Clipped from MACROMS 2009 simulation of the Middle Atlantic Bight (implemented by Wen Long).**

[baymouth_grid_macroms.nc](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/baymouth_grid_macroms.nc)

**ROMS predictions (for LTRANS v.2b input). Files are ~777 MB each. Clipped from MACROMS 2009 simulation of the Middle Atlantic Bight (implemented by Wen Long).**

[clipped_macroms_his_0003.nc](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/clipped_macroms_his_0003.nc)  
[clipped_macroms_his_0004.nc](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/clipped_macroms_his_0004.nc)  
[clipped_macroms_his_0005.nc](http://northweb.hpl.umces.edu/LTRANS/LTRANS-v2/clipped_macroms_his_0005.nc) 	
