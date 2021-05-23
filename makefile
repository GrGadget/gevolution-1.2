# programming environment
COMPILER     := c++
INCLUDE      := $(shell pkg-config --cflags latfield fftw3 hdf5 gsl ompi)
LIB          := $(shell pkg-config --libs latfield fftw3 hdf5 gsl ompi)
HPXCXXLIB    := $(shell pkg-config --libs --cflags healpix)

# Boost
INCLUDE      += -I$(BOOST_ROOT)/include
LIB          += -L$(BOOST_ROOT)/lib -lboost_mpi -lboost_serialization

# target and source
EXEC         := gevolution
SOURCE       := main.cpp
HEADERS      := $(wildcard *.hpp)
SOURCES      := main.cpp lccat.cpp lcmap.cpp hibernation.cpp velocity.cpp \
    background.cpp tools.cpp output.cpp gevolution.cpp radiation.cpp \
    parser.cpp ic_basic.cpp ic_prevolution.cpp ic_read.cpp debugger.cpp \
    newtonian_pm.cpp Particles_gevolution.cpp
VERSION      := version.h
OBJS         := hibernation.o main.o velocity.o background.o tools.o \
    output.o gevolution.o radiation.o parser.o ic_basic.o ic_prevolution.o \
    ic_read.o debugger.o Particles_gevolution.o

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
OPT          := -g -std=c++17

$(EXEC): $(OBJS) $(HEADERS) makefile $(VERSION)
	$(COMPILER) $(OBJS) -o $@ $(OPT) $(DLATFIELD2) $(DGEVOLUTION) $(INCLUDE) $(LIB)

$(OBJS) : %.o : %.cpp $(VERSION)
	$(COMPILER) -c $^ $(INCLUDE) $(DGEVOLUTION) $(OPT)

lccat: lccat.cpp
	$(COMPILER) $< -o $@ $(OPT) $(DGEVOLUTION) $(INCLUDE)
	
lcmap: lcmap.cpp
	$(COMPILER) $< -o $@ $(OPT) -fopenmp $(DGEVOLUTION) $(INCLUDE) $(LIB) $(HPXCXXLIB)

$(VERSION) : git_version.sh version.template
	bash $< . . version.template

clean:
	-rm -f $(EXEC) lccat lcmap $(VERSION) $(OBJS)

format:
	-clang-format -i $(HEADERS) $(SOURCES)

.PHONY: clean format
