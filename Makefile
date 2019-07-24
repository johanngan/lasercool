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
$(srcdir)/fundconst:\
$(srcdir)/optmol:\
$(srcdir)/swapcool

# swapmotion needs extra linking for GSL integration
PROG = optical_molasses swapint
BINS = $(addprefix $(bindir)/, $(PROG))
ARCHIVES = libreadcfg.a libiotag.a libfundconst.a
LIBS = $(addprefix $(libdir)/, $(ARCHIVES))

.PHONY: all clean libs readcfg iotag fundconst optmol swapint swapmotion swapcool
all: $(LIBS) $(BINS) swapmotion
libs: $(LIBS)
readcfg: $(libdir)/libreadcfg.a
iotag: $(libdir)/libiotag.a
fundconst: $(libdir)/libfundconst.a
optmol: $(bindir)/optical_molasses
swapint: $(bindir)/swapint
swapmotion: $(bindir)/swapmotion
swapcool: swapint swapmotion

$(BINS):
	$(LD) $(ALL_LFLAGS) $^ -L$(libdir) -lreadcfg -liotag -lfundconst -o $@

$(bindir)/swapmotion:
	$(LD) $(ALL_LFLAGS) $^ -L$(libdir) -lreadcfg -liotag -lgsl -lgslcblas -lm -o $@

$(bindir)/optical_molasses: \
$(builddir)/optical_molasses.o \
$(builddir)/constants.o \
$(builddir)/PhysicalParams.o \
$(libdir)/libreadcfg.a \
$(libdir)/libiotag.a \
$(libdir)/libfundconst.a

$(bindir)/swapint: \
$(builddir)/swapint.o \
$(builddir)/HInt.o \
$(builddir)/HSwap.o \
$(libdir)/libreadcfg.a \
$(libdir)/libiotag.a \
$(libdir)/libfundconst.a

$(bindir)/swapmotion: \
$(builddir)/swapmotion.o \
$(builddir)/HMotion.o \
$(builddir)/HSwap.o \
$(builddir)/DensMatHandler.o \
$(libdir)/libreadcfg.a \
$(libdir)/libiotag.a \
$(libdir)/libfundconst.a

$(builddir)/optical_molasses.o: optical_molasses.cpp mathutil.hpp RandProcesses.hpp
$(builddir)/PhysicalParams.o: PhysicalParams.cpp mathutil.hpp
$(builddir)/swapint.o: swapint.cpp timestepping.hpp
$(builddir)/swapmotion.o: swapmotion.cpp timestepping.hpp

$(builddir)/optical_molasses.o:
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) -I$(vendordir)/pcg-cpp-0.98/include $< -o $@

$(builddir)/PhysicalParams.o \
$(builddir)/swapint.o \
$(builddir)/swapmotion.o:
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) $< -o $@

$(builddir)/%.o: %.cpp
	$(CC) -c $(ALL_CFLAGS) -I$(includedir) $< -o $@

$(libdir)/lib%.a: $(builddir)/%.o
	$(AR) $(ARFLAGS) $@ $<

clean:
	rm -f $(bindir)/* $(builddir)/*.o $(libdir)/*.a
