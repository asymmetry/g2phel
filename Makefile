OS   := $(shell uname)
ARCH := $(shell arch)

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler))

LIBDIR       := $(ANALYZER)
LIBET        := $(CODA)/Linux/lib/libet.a
HALLALIBS    := -L$(LIBDIR) -lHallA -ldc -lscaler

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

INCLUDES      = $(ROOTCFLAGS) $(addprefix -I, $(INCDIRS))

CXX           = g++
CXXFLAGS      = -O3
LD            = g++
LDFLAGS       = -O3 
DEFINES       = -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC -fno-strict-aliasing

CXXFLAGS     += $(DEFINES) $(INCLUDES)
LIBS         += $(ROOTLIBS) $(HALLALIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(SYSLIBS)

PROGRAMS = extract ring tir align gen_hel

SRCDIR := src

all: $(PROGRAMS)

extract: $(SRCDIR)/extract.o
	$(LD) -g $(CXXFLAGS) -o $@ $< $(LIBS)

ring: $(SRCDIR)/ring.o
	$(LD) -g $(CXXFLAGS) -o $@ $< $(ROOTLIBS)

tir: $(SRCDIR)/tir.o
	$(LD) -g $(CXXFLAGS) -o $@ $< $(ROOTLIBS)

align: $(SRCDIR)/align.o
	$(LD) -g $(CXXFLAGS) -o $@ $< $(ROOTLIBS)

gen_hel: $(SRCDIR)/gen.o
	$(LD) -g $(CXXFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(SRCDIR)/*.o $(PROGRAMS)

release:
	tar cvzf helicity.tar.gz HISTORY Makefile src --exclude=*.o

$(SRCDIR)/%.o:	$(SRCDIR)/%.cxx
	$(CXX) $(CXXFLAGS) -o $@ -c $<
