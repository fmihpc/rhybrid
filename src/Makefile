USE_NODE_UE := true
USE_EDGE_J := true
USE_B_INITIAL := true
USE_B_CONSTANT := false
USE_RESISTIVITY := true
USE_BACKGROUND_CHARGE_DENSITY := false
USE_OUTER_BOUNDARY_ZONE := false
USE_GRID_CONSTRAINT_COUNTERS := false
WRITE_GRID_TEMPORAL_AVERAGES := true
USE_DETECTORS := false
USE_SHOCKTUBE_TEST_CONFIGURATION := false
USE_NEW_VARIBLE_HANDLING := true

include ../../../Makefile.${ARCH}

# Names of files that should be compiled and inserted to lib${SIM}.a:
OBJS = detectors.o hybrid.o hybrid_propagator.o operator_userdata.o particle_accumulator.o particle_injector.o particle_species.o register_objects.o user.o

# Set RHybrid compilation options/macros

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

ifeq ($(USE_RESISTIVITY),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_RESISTIVITY
endif

ifeq ($(USE_BACKGROUND_CHARGE_DENSITY),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_BACKGROUND_CHARGE_DENSITY
endif

ifeq ($(USE_OUTER_BOUNDARY_ZONE),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_OUTER_BOUNDARY_ZONE
endif

ifeq ($(USE_GRID_CONSTRAINT_COUNTERS),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_GRID_CONSTRAINT_COUNTERS
endif

ifeq ($(WRITE_GRID_TEMPORAL_AVERAGES),true)
CXXFLAGS := $(CXXFLAGS) -DWRITE_GRID_TEMPORAL_AVERAGES
endif

ifeq ($(USE_DETECTORS),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_DETECTORS
endif

ifeq ($(USE_SHOCKTUBE_TEST_CONFIGURATION),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_SHOCKTUBE_TEST_CONFIGURATION
endif

ifeq ($(USE_NEW_VARIBLE_HANDLING),true)
CXXFLAGS := $(CXXFLAGS) -DUSE_NEW_VARIBLE_HANDLING
endif

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

# Compile information passed to user.cpp
COMPILEINFO := \""DATE: $(shell date --iso-8601=seconds)\nUSER: $(shell whoami)\nHOST: $(shell hostname)\nDIR: $(shell pwd)\nCMP: $(shell echo ${CMP})\nCXXFLAGS: $(shell echo ${CXXFLAGS})\nMATHFLAGS: $(shell echo ${MATHFLAGS})\nFLAGS: $(shell echo ${FLAGS})\nINCS: $(shell echo ${INCS})\nINCS_REG: $(shell echo ${INCS_REG})\""

# Compilation rules

.cpp.o:
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c $< ${INCS}

detectors.o: detectors.cpp detectors.h hybrid.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c detectors.cpp ${INCS}

hybrid.o: hybrid.cpp hybrid.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c hybrid.cpp ${INCS}

hybrid_propagator.o: hybrid_propagator.cpp hybrid.h hybrid_propagator.h  particle_definition.h magnetic_field.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c hybrid_propagator.cpp ${INCS}

operator_userdata.o: operator_userdata.cpp hybrid.h operator_userdata.h particle_species.h hybrid_propagator.h particle_definition.h magnetic_field.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c operator_userdata.cpp ${INCS_REG}

particle_accumulator.o: particle_accumulator.cpp hybrid.h particle_accumulator.h particle_definition.h particle_species.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_accumulator.cpp ${INCS}

particle_injector.o: particle_injector.cpp particle_injector.h particle_definition.h particle_species.h hybrid.h neutral_profiles.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_injector.cpp ${INCS_REG}

particle_species.o: particle_species.cpp particle_species.h hybrid.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c particle_species.cpp ${INCS_REG}

register_objects.o: register_objects.cpp particle_boundary_cond_hybrid.h particle_definition.h particle_species.h operator_userdata.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -c register_objects.cpp ${INCS_REG}

user.o: user.cpp hybrid.h hybrid_propagator.h particle_definition.h particle_species.h particle_accumulator.h particle_propagator_boris_buneman.h particle_injector.h particle_list_hybrid.h operator_userdata.h diagnostics.h detectors.h magnetic_field.h
	${CMP} ${CXXFLAGS} ${MATHFLAGS} ${FLAGS} -DCOMPILEINFO=${COMPILEINFO} -c user.cpp ${INCS} ${INCS_REG}
