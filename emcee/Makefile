SRCDIR = ./src
OBJDIR = ./build
CXXFLAGS = -I${SRCDIR} -std=c++11 -O3 
CPP_FILES := $(wildcard $(SRCDIR)/*emcee*.cc)
OBJ_FILES := $(addprefix ${OBJDIR}/,$(notdir $(CPP_FILES:.cc=.o)))
TARGETS = example_sinc example_lighthouse example_gaussian example_rosenbrock

all: ${TARGETS}


build/%.o: ${SRCDIR}/%.cc
	@mkdir -p $(@D)
	g++ $(CXXFLAGS) -c -o $@ $<

example_lighthouse: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^ $(SRCDIR)/$@.cc

example_sinc: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^ $(SRCDIR)/$@.cc

example_gaussian: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^ $(SRCDIR)/$@.cc

example_rosenbrock: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^ $(SRCDIR)/$@.cc

clean:
	rm -f ${TARGETS}
	rm -rf ${OBJDIR}

