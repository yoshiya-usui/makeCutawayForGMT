CXX           = icpc
CC            = icpc
COMP_ROOT     = /opt/intel/composer_xe_2015.5.223
COMP_LIB      = ${COMP_ROOT}/lib/intel64
CXXFLAGS      = -O3 \
                -openmp \
                -D_LINUX \
                -DNDEBUG \
                -I${COMP_ROOT}/include
DEST          = ./
LDFLAGS       = -L${COMP_LIB}
LIBS          = 
OBJS          = main.o \
                MeshData.o \
                MeshDataTetraElement.o \
                MeshDataBrickElement.o \
                ResistivityBlock.o \
                Util.o
PROGRAM       = makeCutawayForGMT

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) $(LDFLAGS) $(LIBS) -o $(PROGRAM)

clean:;		rm -f *.o *~ $(PROGRAM)
