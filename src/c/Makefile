SRC = reg_xyz.cpp
HEADERS = reg_xyz.h
CXX = g++

system := $(shell uname)
ifeq ($(system),Linux)
CFLAGS=-fPIC
LIBS=-lgsl -lm -lgslcblas
else
INC=-I/sw/include
LIBDIR=-L/sw/lib
LIBS=-lgsl -lm
endif

all: reg_xyz.so

reg_xyz.so: $(SRC) $(HEADERS)
	PKG_LIBS="$(LIBDIR) $(LIBS)" R CMD SHLIB -o reg_xyz.so $(SRC)

clean:
	-rm *.o reg_xyz.so

# END
