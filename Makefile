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
CXXFLAGS      = -O3 -Wall -Woverloaded-virtual -fPIC -fno-strict-aliasing
CXXFLAGS     += $(INCLUDES)

PROGRAMS = hel_decode hel_ring hel_happex hel_tir hel_alignr hel_alignh hel_insert

SRCDIR := src

all: $(PROGRAMS)

hel_decode: $(SRCDIR)/decode.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS) $(HALLALIBS) $(CONFLIBS)

hel_ring: $(SRCDIR)/ring.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

hel_happex: $(SRCDIR)/ring.cxx
	$(CXX) $(CXXFLAGS) -DHAPPEX -o $@ $< $(ROOTLIBS) $(CONFLIBS)

hel_tir: $(SRCDIR)/tir.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

hel_alignr: $(SRCDIR)/align.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS) $(CONFLIBS)

hel_alignh: $(SRCDIR)/align.cxx
	$(CXX) $(CXXFLAGS) -DHAPPEX -o $@ $< $(ROOTLIBS) $(CONFLIBS)

hel_insert: $(SRCDIR)/insert.cxx
	$(CXX) $(CXXFLAGS) -o $@ $< $(ROOTLIBS) $(HALLALIBS) $(CONFLIBS)

clean:
	rm -f $(PROGRAMS)

release:
	tar cvzf helicity.tar.gz README Makefile src helicity.sh --exclude=*.o
