USE_NODE_UE := true
USE_EDGE_J := true
USE_B_INITIAL := true
USE_B_CONSTANT := false
USE_XMIN_BOUNDARY := false
USE_CONIC_INNER_BOUNDARY := false
USE_RESISTIVITY := true
USE_ECUT := true
USE_MAXVW := true
WRITE_POPULATION_AVERAGES := true
USE_DETECTORS := false

include ../../../Makefile.${ARCH}

# Names of files that should be compiled and 
# inserted to lib${SIM}.a:

OBJS = register_objects.o user.o hybrid_propagator.o hybrid.o\
	particle_accumulator.o particle_injector.o\
	operator_userdata.o particle_species.o detectors.o

ifeq ($(USE_NODE_UE),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_NODE_UE
endif

ifeq ($(USE_EDGE_J),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_EDGE_J
endif

ifeq ($(USE_B_INITIAL),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_B_INITIAL
endif

ifeq ($(USE_B_CONSTANT),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_B_CONSTANT
endif

ifeq ($(USE_XMIN_BOUNDARY),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_XMIN_BOUNDARY
endif

ifeq ($(USE_CONIC_INNER_BOUNDARY),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_CONIC_INNER_BOUNDARY
endif

ifeq ($(USE_RESISTIVITY),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_RESISTIVITY
endif

ifeq ($(USE_ECUT),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_ECUT
endif

ifeq ($(USE_MAXVW),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_MAXVW
endif

ifeq ($(WRITE_POPULATION_AVERAGES),true)
CXXFLAGS := $(CXXFLAGS) -DWRITE_POPULATION_AVERAGES
endif

ifeq ($(USE_DETECTORS),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_DETECTORS
endif

#override CXXFLAGS += -std=gnu++0x

# Compile information (date___user___host___folder___cxxflags)
COMPILEINFO := \"$(shell date --iso-8601=seconds)___$(shell whoami)___$(shell hostname)___$(shell pwd)___$(shell echo ${CXXFLAGS})\"
# Remove spaces and add newlines
COMPILEINFO := \"$(shell echo $(COMPILEINFO) | sed 's/ //g' | sed 's/___/\\n/g' | sed 's/-D/\\n/g')\"

# Make targets

default:
	${MAKE} lib${SIM}.a

clean:
	rm -rf *.o *.a *~
	rm -f ../lib/lib${SIM}.a

lib${SIM}.a: ${OBJS}
	${AR} r lib${SIM}.a ${OBJS}
	ln -f -s ${CURDIR}/lib${SIM}.a -t ../../lib

# Dependencies

override INCS+=${INC} ${INC_BOOST} ${INC_PARGRID} ${INC_VLSV} ${INC_ZOLTAN}
override INCS+=-I../../include

INCS_REG=${INCS} -I../../particleinjector -I../../dataoperator
INCS_REG+=-I../../particlepropagator -I../../gridbuilder

DEPS_ACCUM=particle_definition.h particle_species.h hybrid.h particle_accumulator.h particle_accumulator.cpp
DEPS_REG_OBJS=register_objects.cpp
DEPS_SPECIES=particle_species.h particle_species.cpp
DEPS_ADV_PROP=hybrid.h hybrid_propagator.h hybrid_propagator.cpp
DEPS_EX_ADV=hybrid.h hybrid.cpp
DEPS_OP_USER=operator_userdata.h operator_userdata.cpp
DEPS_INJECTOR=particle_definition.h particle_species.h particle_injector.h particle_injector.cpp
DEPS_USER=${DEPS_ACCUM} ${DEPS_SPECIES} ${DEPS_EX_ADV} ${DEPS_INJECTOR} particle_propagator_boris_buneman.h ../../include/user.h user.cpp particle_list_hybrid.h

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INCS}

register_objects.o: ${DEPS_REG_OBJS}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c register_objects.cpp ${INCS_REG}

user.o: ${DEPS_USER}
	${CMP} ${CXXFLAGS} -DCOMPILEINFO=${COMPILEINFO} ${MATHFLAGS} ${FLAGS} -c user.cpp ${INCS} ${INCS_REG}

hybrid_propagator.o: ${DEPS_ADV_PROP}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c hybrid_propagator.cpp ${INCS}

hybrid.o: ${DEPS_EX_ADV}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c hybrid.cpp ${INCS}

operator_userdata.o: ${DEPS_OP_USER}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c operator_userdata.cpp ${INCS_REG}

particle_accumulator.o: ${DEPS_ACCUM}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_accumulator.cpp ${INCS}

particle_injector.o: ${DEPS_INJECTOR}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_injector.cpp ${INCS_REG}	

particle_species.o: ${DEPS_SPECIES}
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_species.cpp ${INCS_REG}

detectors.o: detectors.h detectors.cpp
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c detectors.cpp ${INCS}
