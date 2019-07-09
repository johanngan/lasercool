SHELL = /bin/sh
AR = ar
ARFLAGS = rvs
CC = g++
CFLAGS =
ALL_CFLAGS = -std=c++14 -O3 -flto -Wall -Wextra $(CFLAGS)
LD = g++
LFLAGS =
ALL_LFLAGS = -O3 -flto $(LFLAGS)

prefix = .
bindir = $(prefix)/bin
builddir = $(prefix)/build
includedir = $(prefix)/include
libdir = $(prefix)/lib
srcdir = $(prefix)/src
vendordir = $(prefix)/vendor

VPATH = \
$(includedir)/lasercool:\
$(srcdir)/iotag:\
$(srcdir)/readcfg:\
$(srcdir)/optmol:\
$(srcdir)/swapcool

# swapmotion needs extra linking for GSL integration
PROG = optical_molasses swapint
BINS = $(addprefix $(bindir)/, $(PROG))
ARCHIVES = libreadcfg.a libiotag.a
LIBS = $(addprefix $(libdir)/, $(ARCHIVES))

.PHONY: all clean libs readcfg iotag optmol swapint swapmotion swapcool
all: $(LIBS) $(BINS) swapmotion
libs: $(LIBS)
readcfg: $(libdir)/libreadcfg.a
iotag: $(libdir)/libiotag.a
optmol: readcfg iotag $(bindir)/optical_molasses
swapint: readcfg iotag $(bindir)/swapint
swapmotion: readcfg iotag $(bindir)/swapmotion
swapcool: swapint swapmotion

$(BINS):
	$(LD) $(ALL_LFLAGS) $^ -L$(libdir) -lreadcfg -liotag -o $@

$(bindir)/swapmotion:
	$(LD) $(ALL_LFLAGS) $^ -L$(libdir) -lreadcfg -liotag -lgsl -lgslcblas -lm -o $@

$(bindir)/optical_molasses: \
$(builddir)/optical_molasses.o \
$(builddir)/constants.o \
$(builddir)/PhysicalParams.o

$(bindir)/swapint: \
$(builddir)/swapint.o \
$(builddir)/HInt.o \
$(builddir)/HSwap.o

$(bindir)/swapmotion: \
$(builddir)/swapmotion.o \
$(builddir)/HMotion.o \
$(builddir)/HSwap.o \
$(builddir)/DensMatHandler.o

$(builddir)/optical_molasses.o: optical_molasses.cpp mathutil.hpp RandProcesses.hpp
$(builddir)/PhysicalParams.o: PhysicalParams.cpp mathutil.hpp
$(builddir)/swapint.o: swapint.cpp timestepping.hpp
$(builddir)/swapmotion.o: swapmotion.cpp timestepping.hpp HMotion.hpp DensMatHandler.hpp
$(builddir)/HMotion.o: HMotion.cpp DensMatHandler.hpp

$(builddir)/optical_molasses.o:
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) -I$(vendordir)/pcg-cpp-0.98/include $< -o $@

$(builddir)/PhysicalParams.o \
$(builddir)/swapint.o \
$(builddir)/swapmotion.o \
$(builddir)/HMotion.o:
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) $< -o $@

$(builddir)/%.o: %.cpp
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) $< -o $@

$(libdir)/lib%.a: $(builddir)/%.o
	$(AR) $(ARFLAGS) $@ $<

clean:
	rm -f $(bindir)/* $(builddir)/*.o $(libdir)/*.a
