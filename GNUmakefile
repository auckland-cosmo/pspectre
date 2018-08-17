# SpectRE - A Spectral Code for Reheating
# Copyright (C) 2009-2010 Hal Finkel
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES,
# INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
# AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
# THE AUTHORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

ifeq ($(USE_ICC),yes)
GXX=icc
GXX_OPT_FLAGS = $(shell sh ./testcompilecull.sh $(GXX) -openmp -ipo -xHost)
GXX_WARN_FLAGS =
GXX_EXTRA_LIBS = -lstdc++
else
GXX=g++
ifneq ($(wildcard $(shell which g++-4 2>/dev/null | grep '^/')),)
GXX=g++-4
endif

ifneq ($(NO_M64),yes)
M64_FLAG = -m64
else
M64_FLAG =
endif

GXX_OPT_FLAGS = $(shell sh ./testcompilecull.sh $(GXX) -fopenmp -ftree-vectorize $(M64_FLAG) -march=native)
GXX_WARN_FLAGS = -Wall
GXX_EXTRA_LIBS =
endif

ifneq ($(USE_MKL),yes)
ifneq ($(wildcard $(shell which fftwl-wisdom 2>/dev/null | grep '^/')),)
USE_LD = yes
endif
endif

ifneq ($(wildcard private/*.cpp),)
ifneq ($(NO_PRIV),yes)
CPPFLAGS_PRIV=-DHAVE_PRIVATE
PRIV_SRCS=$(wildcard private/*.cpp)
PRIV_CLEAN=private/*.o
endif
endif

SRCS= \
	field.cpp \
	model.cpp \
	integrator.cpp \
	verlet.cpp \
	rk4.cpp \
	v_integrator.cpp \
	nonlinear_transformer.cpp \
	le_style_initializer.cpp \
	defrost_style_initializer.cpp \
	slice_output_manager.cpp \
	slice_outputter.cpp \
	spectra_outputter.cpp \
	twoptcorr_outputter.cpp \
	stats_outputter.cpp \
	energy_outputter.cpp \
	grid_funcs.cpp \
	grad_computer.cpp \
	gpot_computer.cpp \
	pspectre.cpp \
	$(PRIV_SRCS)

rel: pspectre
profile: pspectre-pg
debug: pspectre-dbg

ifeq ($(USE_LD),yes)
CPPFLAGS_EXTRA = -DUSE_LD
endif

ifeq ($(USE_MKL),yes)
CPPFLAGS_EXTRA = -DUSE_MKL

ifneq ($(wildcard $(MKLROOT)/include/fftw/fftw3_mkl.h),)
CPPFLAGS_EXTRA += -DHAS_FFTW3_MKL_H
MKL_IS_11 = yes
else
CPPFLAGS_EXTRA += -DMKL_NO_DCT
MKL_IS_11 = no
endif

MKL_FFTW_LIB_ARCH = 32
ifneq ($(shell $(GXX) -dumpmachine | sed -n '/^x86_64-/p;/^amd64-/p'),)
MKL_FFTW_LIB_ARCH = em64t
endif
ifneq ($(shell $(GXX) -dumpmachine | sed -n '/^ia64-/p'),)
MKL_FFTW_LIB_ARCH = 64
endif

libfftw3xc_intel.a:
	mkdir -p fftw3xc-build/lib fftw3xc-build/lib/$(MKL_FFTW_LIB_ARCH) fftw3xc-build/interfaces
	ln -s $(MKLROOT)/include fftw3xc-build/include
	cp -r $(MKLROOT)/interfaces/fftw3xc fftw3xc-build/interfaces/
	$(MAKE) -C fftw3xc-build/interfaces/fftw3xc lib$(MKL_FFTW_LIB_ARCH)
	find fftw3xc-build -name $@ -exec cp {} $@ \;
	rm -rf fftw3xc-build

libfftw3xc_gnu.a:
	mkdir -p fftw3xc-build/lib fftw3xc-build/lib/$(MKL_FFTW_LIB_ARCH) fftw3xc-build/interfaces
	ln -s $(MKLROOT)/include fftw3xc-build/include
	cp -r $(MKLROOT)/interfaces/fftw3xc fftw3xc-build/interfaces/
	$(MAKE) -C fftw3xc-build/interfaces/fftw3xc lib$(MKL_FFTW_LIB_ARCH) compiler=gnu
	find fftw3xc-build -name $@ -exec cp {} $@ \;
	rm -rf fftw3xc-build

ifeq ($(USE_ICC),yes)
MKL_FFTW_WRAP_LIB = libfftw3xc_intel.a
else
MKL_FFTW_WRAP_LIB = libfftw3xc_gnu.a
endif
endif

CPPFLAGS_EXTRA += $(CPPFLAGS_PRIV)

ifneq ($(shell sh ./testcompilecull.sh $(GXX) -fmudflap -lmudflap),)
HAVE_MUDFLAP=yes
debug-mudflap: pspectre-dbg-mf
else
HAVE_MUDFLAP=no
endif

GXXFLAGS = -O3 $(GXX_OPT_FLAGS) $(CPPFLAGS_EXTRA) $(GXX_WARN_FLAGS)
GXXFLAGS_PG = -O3 $(GXX_OPT_FLAGS) $(CPPFLAGS_EXTRA) $(GXX_WARN_FLAGS) -pg
GXXFLAGS_DBG = -O0 $(CPPFLAGS_EXTRA) $(GXX_WARN_FLAGS) -g
GXXFLAGS_DBG_MF = -O0 $(CPPFLAGS_EXTRA) $(GXX_WARN_FLAGS) -g -fmudflap

DOXYGEN_OUTPUT = documentation/html/index.html documentation/latex/refman.tex

ifneq ($(NO_DEP_GEN),yes)
mk/%.d: %.cpp
	@echo performing automated dependency generation for $<
	@set -e; mkdir -p $(dir $@); rm -f $@; \
	$(GXX) -M $(GXXFLAGS) $< 2>/dev/null > $@.$$$$; \
	sed 's,\($*\)\.o[ :]*,\1.o $@ $(DOXYGEN_OUTPUT) : ,g' < $@.$$$$ > $@; \
	echo >> $@; \
	sed 's,\($*\)\.o[ :]*,\1-pg.o $@ $(DOXYGEN_OUTPUT) : ,g' < $@.$$$$ >> $@; \
	echo >> $@; \
	sed 's,\($*\)\.o[ :]*,\1-dbg.o $@ $(DOXYGEN_OUTPUT) : ,g' < $@.$$$$ >> $@; \
	echo >> $@; \
	sed 's,\($*\)\.o[ :]*,\1-dbg-mf.o $@ $(DOXYGEN_OUTPUT) : ,g' < $@.$$$$ >> $@; \
	rm -f $@.$$$$
endif

-include $(addprefix mk/,$(SRCS:.cpp=.d))

%.o: %.cpp
	$(GXX) -c $(GXXFLAGS) -o $@ $<

%-pg.o: %.cpp
	$(GXX) -c $(GXXFLAGS_PG) -o $@ $<

%-dbg.o: %.cpp
	$(GXX) -c $(GXXFLAGS_DBG) -o $@ $<

ifeq ($(HAVE_MUDFLAP),yes)
%-dbg-mf.o: %.cpp
	$(GXX) -c $(GXXFLAGS_DBG_MF) -o $@ $<
endif

ifeq ($(USE_MKL),yes)
FFT_LIBS = $(MKL_FFTW_WRAP_LIB)

ifeq ($(USE_ICC),yes)
ifeq ($(MKL_IS_11),yes)
ifeq ($(MKL_FFTW_LIB_ARCH),32)
FFT_LIBS += -lmkl_intel
else
FFT_LIBS += -lmkl_intel_lp64
endif
else
FFT_LIBS += -lmkl_intel_ilp64
endif
FFT_LIBS += -lmkl_intel_thread
else
ifeq ($(MKL_IS_11),yes)
FFT_LIBS += -lmkl_gf
else
FFT_LIBS += -lmkl_gf_ilp64
endif
FFT_LIBS += -lmkl_gnu_thread
endif

FFT_LIBS += -lmkl_core

ifneq ($(MKL_IS_11),yes)
FFT_LIBS += -lguide
endif

else
FFT_LIBS = -lfftw3 -lfftw3_threads

ifeq ($(USE_LD),yes)
FFT_LIBS += -lfftw3l -lfftw3l_threads
endif
endif

pspectre: $(SRCS:.cpp=.o) $(MKL_FFTW_WRAP_LIB)
	$(GXX) $(GXXFLAGS) -o $@ $^ $(FFT_LIBS) $(GXX_EXTRA_LIBS)

pspectre-pg: $(SRCS:.cpp=-pg.o) $(MKL_FFTW_WRAP_LIB)
	$(GXX) $(GXXFLAGS_PG) -o $@ $^ $(FFT_LIBS) $(GXX_EXTRA_LIBS)

pspectre-dbg: $(SRCS:.cpp=-dbg.o) $(MKL_FFTW_WRAP_LIB)
	$(GXX) $(GXXFLAGS_DBG) -o $@ $^ $(FFT_LIBS) $(GXX_EXTRA_LIBS)

ifeq ($(HAVE_MUDFLAP),yes)
pspectre-dbg-mf: $(SRCS:.cpp=-dbg-mf.o) $(MKL_FFTW_WRAP_LIB)
	$(GXX) $(GXXFLAGS_DBG_MF) -o $@ $^ $(FFT_LIBS) $(GXX_EXTRA_LIBS) -lmudflap
endif

doc: documentation/html/index.html documentation/latex/refman.pdf

$(DOXYGEN_OUTPUT): Doxyfile $(SRCS)
	doxygen

documentation/latex/refman.pdf: documentation/latex/refman.tex
	$(MAKE) -C documentation/latex pdf
	mv documentation/latex/refman.pdf documentation/latex/refman.pdf.save
	$(MAKE) -C documentation/latex clean
	mv documentation/latex/refman.pdf.save documentation/latex/refman.pdf

clean:
	rm -f *.o $(PRIV_CLEAN) *.a pspectre pspectre-pg pspectre-dbg pspectre-dbg-mf
	rm -rf mk fftw3xc-build

clean-doc:
	rm -rf documentation

clean-all: clean clean-doc

dist:
	$(MAKE) clean-all
	$(MAKE) doc
	find mk >> ../pspectre-dist-exclude
	if test -d .svn; then find .svn >> ../pspectre-dist-exclude; fi
	if test -d private; then find private >> ../pspectre-dist-exclude; fi
	for odir in output-*; do if test -d "$$odir"; then find $$odir >> ../pspectre-dist-exclude; fi; done
	cd .. && tar -X pspectre-dist-exclude -czvf pspectre-$(shell date '+%Y%m%d').tar.gz pspectre
	rm -f ../pspectre-dist-exclude

