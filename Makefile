#-----------------------------------------------------------------------
#   Configuration of compile
#-----------------------------------------------------------------------
# Variables that must be specified by users
# (i) Variables related to FDPS
FDPS_LOC = ../fdps/FDPS/
FDPS_INC = -I$(FDPS_LOC)/src 
FDPS_INC += -I$(FDPS_LOC)/src/c_interface/headers
FDPS_C_IF_GENERATOR = $(FDPS_LOC)/scripts/gen_c_if.py


# (ii) Variables to specify compilers and compile options
# Serial or OpenMP cases
FC=gfortran
CXX=g++
# MPI case
#FC=mpif90
#CXX=mpic++
# [Option 1] w/o optimization
#FCFLAGS = -std=f2003 -O0 -Wall
CXXFLAGS = -Wall -Wextra -ftrapv -fexceptions -g3 $(FDPS_INC)
# [Option 2] w/ optimization 
FCFLAGS = -std=f2003 -O3 -ffast-math -funroll-loops -finline-functions
#CXXFLAGS = -O3 -ffast-math -funroll-loops $(FDPS_INC)
CXXFLAGS = -g -funroll-loops $(FDPS_INC)
LDFLAGS = -lgfortran 
# OpenMP options
#FCFLAGS  += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
#CXXFLAGS += -DPARTICLE_SIMULATOR_THREAD_PARALLEL -fopenmp
# MPI options
#FCFLAGS  += -DPARTICLE_SIMULATOR_MPI_PARALLEL
#CXXFLAGS += -DPARTICLE_SIMULATOR_MPI_PARALLEL

# fdps-autotest-set-vars (DO NOT CHANGE THIS LINE)

#-----------------------------------------------------------------------
#   Source files
#-----------------------------------------------------------------------

%.o : %.c
	$(CC) $(CFLAGS) -c $<
%.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $<

HDR_USER_DEFINED_TYPE = user_defined.h
SRC_USER = user_defined.c \
	   c_main.c 
SRC_CXX = FDPS_ftn_if.cpp \
	  FDPS_Manipulators.cpp \
	  main.cpp

OBJ_USER = $(SRC_USER:c=o)
OBJ_CXX	 = $(SRC_CXX:cpp=o)
OBJ	 = $(OBJ_USER) $(OBJ_CXX)


FDPS_c_if.h $(SRC_CXX): $(HDR_USER_DEFINED_TYPE) Makefile
	$(FDPS_C_IF_GENERATOR) user_defined.h --output ./
FDPS_cr_if.cr:   FDPS_ftn_if.cpp  convert_f90_if_to_crystal.rb
	ruby convert_f90_if_to_crystal.rb FDPS_ftn_if.cpp > FDPS_cr_if.cr

$(OBJ_USER): FDPS_c_if.h 

user_defined.h: user_defined.cr
	ruby convert_crystal_struct_to_c.rb user_defined.cr >user_defined.h
libcrmainx.so: crmain.cr user_defined.cr FDPS_vector.cr FDPS_cr_if.cr FDPS_types.cr Makefile
	crystal build  crmain.cr  --threads 1 --single-module --link-flags="-shared" -o libcrmainx.so

CROBJS =   FDPS_ftn_if.o FDPS_Manipulators.o crmain.o
CRLIBS = libcrmainx.so
CLIBS = -levent -ldl -lpcre 
fdpscr:  $(CROBJS) $(CRLIBS) Makefile
	$(CXX) $(CXXFLAGS) $(CROBJS)    -o fdpscr -L. -lcrmainx  $(CLIBS)
crmain.o: crmain.cpp FDPS_Manipulators.h user_defined.h
clean:
	rm -f *.o *.s $(TARGET) *.dat

distclean: clean
	rm -f $(SRC_CXX) FDPS_c_if.h FDPS_Manipulators.h user_defined.hpp 
	rm -rf result


EXPORTDIR = export
EXPORTSRCS = Makefile   convert_crystal_struct_to_f90.rb\
             convert_crystal_struct_to_c.rb crmain.cr\
             FDPS_vector.cr README.md  convert_f90_if_to_crystal.rb\
             user_defined.cr LICENSE	convert_f90_struct_to_crystal.rb\
             shard.yml



exports : $(EXPORTSRCS)
	cp -p $(EXPORTSRCS) $(EXPORTDIR)

git :
	git add --all
	git commit
	git push origin master
