OS   := $(shell uname)
ARCH := $(shell arch)

ifndef ANALYZER
  $(error $$ANALYZER environment variable not defined)
endif

INCDIRS = $(wildcard $(addprefix $(ANALYZER)/, include src hana_decode hana_scaler))

LIBDIR       := $(ANALYZER)
HALLALIBS    := -L$(LIBDIR) -lHallA -ldc -lscaler

ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell root-config --glibs)

INCLUDES      = $(ROOTCFLAGS) $(addprefix -I, $(INCDIRS))

ifndef LIBCONFIG
  $(error $$LIBCONFIG environment variable not defined)
endif
CONFLIBS     := -L$(LIBCONFIG)/lib -lconfig
INCLUDES     += -I$(LIBCONFIG)/include

CXX           = g++
CXXFLAGS      = -O3
LD            = g++
LDFLAGS       = -O3 
DEFINES       = -DLINUXVERS -DHAS_SSTREAM
CXXFLAGS     += -Wall -Woverloaded-virtual -fPIC -fno-strict-aliasing

CXXFLAGS     += $(DEFINES) $(INCLUDES)
LIBS         += $(ROOTLIBS) $(HALLALIBS) $(CONFLIBS) $(SYSLIBS)
GLIBS        += $(ROOTGLIBS) $(HALLALIBS) $(CONFLIBS) $(SYSLIBS)

PROGRAMS = decode ring tir align gen_hel

SRCDIR := src

all: $(PROGRAMS)

decode: $(SRCDIR)/decode.o
	$(LD) -g $(LDFLAGS) -o $@ $< $(LIBS)

ring: $(SRCDIR)/ring.o
	$(LD) -g $(LDFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

tir: $(SRCDIR)/tir.o
	$(LD) -g $(LDFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

align: $(SRCDIR)/align.o
	$(LD) -g $(LDFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

gen_hel: $(SRCDIR)/gen.o
	$(LD) -g $(LDFLAGS) -o $@ $< $(LIBS)

clean:
	rm -f $(SRCDIR)/*.o $(PROGRAMS)

release:
	tar cvzf helicity.tar.gz README Makefile src helicity.sh --exclude=*.o

$(SRCDIR)/%.o:	$(SRCDIR)/%.cxx
	$(CXX) $(CXXFLAGS) -o $@ -c $<
