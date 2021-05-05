# programming environment
COMPILER     := mpic++
INCLUDE      := $(shell pkg-config --cflags latfield fftw3 hdf5 gsl)
LIB          := $(shell pkg-config --libs fftw3 hdf5 gsl)
HPXCXXLIB    := $(shell pkg-config --libs --cflags healpix)

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)
VERSION      := version.h

# mandatory compiler settings (LATfield2)
DLATFIELD2   := -DFFT3D -DHDF5

# optional compiler settings (LATfield2)
#DLATFIELD2   += -DH5_HAVE_PARALLEL
#DLATFIELD2   += -DEXTERNAL_IO # enables I/O server (use with care)
#DLATFIELD2   += -DSINGLE      # switches to single precision, use LIB -lfftw3f

# optional compiler settings (gevolution)
DGEVOLUTION  := -DPHINONLINEAR
DGEVOLUTION  += -DBENCHMARK
DGEVOLUTION  += -DEXACT_OUTPUT_REDSHIFTS
#DGEVOLUTION  += -DVELOCITY      # enables velocity field utilities
#DGEVOLUTION  += -DCOLORTERMINAL
#DGEVOLUTION  += -DCHECK_B
#DGEVOLUTION  += -DHAVE_CLASS    # requires LIB -lclass
#DGEVOLUTION  += -DHAVE_HEALPIX  # requires LIB -lchealpix

# further compiler options
OPT          := -O3 -std=c++11

$(EXEC): $(SOURCE) $(HEADERS) makefile $(VERSION)
	$(COMPILER) $< -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)
	
lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
	
lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

$(VERSION) : git_version.sh version.template
	bash $< . . version.template

clean:
	-rm -f $(EXEC) lccat lcmap $(VERSION)

