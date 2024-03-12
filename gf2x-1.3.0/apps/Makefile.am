# add -mbmi2 to use the pdep_u64 instruction
AM_CXXFLAGS = $(GF2X_CFLAGS)
LDADD = $(GF2X_LIBS) $(GF2X_RPATH_LIBS) -lgmp

bin_PROGRAMS=factor
factor_SOURCES=factor.cpp halfgcd.cpp
noinst_HEADERS=halfgcd.hpp
factor_LDFLAGS=-fopenmp

LOG_COMPILER=$(SHELL) 2>/dev/null
TESTS=check-factor-000.sh check-factor-001.sh	\
      check-factor-002.sh check-factor-003.sh

EXTRA_DIST=$(TESTS)
