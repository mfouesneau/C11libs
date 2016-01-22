SRCDIR = .
OBJDIR = ./build
CXXFLAGS = -I${SRCDIR} -std=c++11 -O3 
CPP_FILES := $(wildcard $(SRCDIR)/*de.cc)
OBJ_FILES := $(addprefix ${OBJDIR}/,$(notdir $(CPP_FILES:.cc=.o)))
TARGETS = de_test


all: ${TARGETS}


build/%.o: ${SRCDIR}/%.cc
	@mkdir -p $(@D)
	g++ $(CXXFLAGS) -c -o $@ $<

de_test: $(OBJ_FILES)
	g++ $(CXXFLAGS) -o $@ $^ $(SRCDIR)/$@.cc

clean:
	rm -rf ${OBJDIR}
	rm -f ${TARGETS}


