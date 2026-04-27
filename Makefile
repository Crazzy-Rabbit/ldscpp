CXX ?= g++
CXXFLAGS ?= -std=c++17 -O2 -Wall -Wextra -pedantic -Iinclude
LDLIBS ?= -lz -lbz2

BIN := bin/ldsc.exe
SRC := src/ldsc.cpp \
       src/common.cpp \
       src/regression.cpp \
       src/sumstats_commands.cpp \
       src/munge.cpp \
       src/annotation.cpp \
       src/ld_score.cpp \
       src/continuous.cpp

.PHONY: all clean

all: $(BIN)

bin:
	if not exist bin mkdir bin

$(BIN): $(SRC) | bin
	$(CXX) $(CXXFLAGS) -o $(BIN) $(SRC) $(LDLIBS)

clean:
	if exist $(BIN) del /Q $(BIN)
