#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
#                                                                       :::
#                          LTRANS Makefile                              :::
#                                                                       :::
#::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::


#==========================================================================
#                         USER-DEFINED OPTIONS                            =
#==========================================================================




#------------------------------------------------
#    Set compiler and flags
#------------------------------------------------
#
#    Turn one of the following on:
IFORT_UMCES := on
GFORTRAN_UMCES :=
GFORTRAN_USGS :=
PGI_USGS   :=

ifdef IFORT_UMCES 
  FC = ifort
  NETCDF_INCDIR = /usr/local/include
  NETCDF_LIBDIR = /usr/local/lib
  FFLAGS = -vec-report0 -fp-model precise -mcmodel=medium -I$(NETCDF_INCDIR)
endif

ifdef GFORTRAN_UMCES
  FC = gfortran
  NETCDF_INCDIR = /data/users/imitchell/local/include
  NETCDF_LIBDIR = /data/users/imitchell/local/lib
  FFLAGS = -march=k8 -ffast-math -fno-cx-limited-range -O3 -funroll-loops --param max-unroll-times=4 -ffree-line-length-none -I$(NETCDF_INCDIR)
endif

ifdef GFORTRAN_USGS
  FC = gfortran
  NETCDF_INCDIR = /usr/local/include
  NETCDF_LIBDIR = /usr/local/lib
  FFLAGS = -Wall  -Warray-bounds -fbacktrace -ffree-line-length-none  -I$(NETCDF_INCDIR)
endif

ifdef PGI_USGS
  FC = pgf90
  NETCDF_INCDIR = /share/apps/netcdf/include
  NETCDF_LIBDIR = /share/apps/netcdf/lib
  FFLAGS := -g -I$(NETCDF_INCDIR)
endif

#------------------------------------------------
#    Set NetCDF Library Locations.
#    If NetCDF was compiled with HDF5, set:
#        HDF5 := on
#    Otherwise, leave blank:
#        HDF5 :=
#------------------------------------------------

HDF5 := on

#==========================================================================
# End of user-defined options. Nothing should be changed below this point =
#==========================================================================

OBJS          = gridcell_module.o interpolation_module.o parameter_module.o \
				point_in_polygon_module.o random_module.o tension_module.o  \
				conversion_module.o hydrodynamic_module.o norm_module.o     \
				boundary_module.o hor_turb_module.o settlement_module.o     \
				ver_turb_module.o behavior_module.o

ifdef HDF5
   ifdef PGI_USGS
        LIBS      = -L$(NETCDF_LIBDIR) -lnetcdf -lnetcdff -L/share/apps/hdf5/lib -lhdf5_hl -lhdf5 -lz -lm -L/share/apps/szip/lib -lsz -lcurl
   else
	LIBS      = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf -lhdf5_hl -lhdf5 -lcurl
   endif
else
	LIBS      = -L$(NETCDF_LIBDIR) -lnetcdff -lnetcdf
endif


LTRANS : $(OBJS)
	@echo "  Compiling LTRANS.f90"
	@$(FC) $(FFLAGS) -o LTRANS.exe LTRANS.f90 $(OBJS) $(LIBS)
	@\rm *.o *.mod
	@echo "  "
	@echo "  Compilation Successfully Completed"
	@echo "  "

%.o: %.f90
	@echo "  Compiling $<"
	@$(FC) $(FFLAGS) -c $<

clean:
	\rm *.o *.mod LTRANS.exe

