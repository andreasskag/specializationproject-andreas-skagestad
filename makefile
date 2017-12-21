
CFLAGS+= -std=c++11 -g -O2

LIBS=-laf
#change this to point to wherever you have arrayfire
AF_PATH = ./../package_portable/arrayfire/#$(HOME)/arrayfire-3/#/usr/local/arrayfire

CUDA_PATH = /usr/local/cuda-9.0/#change this so that it points to your cuda folder
CUDA_LIB_PATH=$(CUDA_PATH)lib64/

INCLUDES=-I $(AF_PATH)include/
LIB_PATH=-L $(AF_PATH)lib/

###CUDA LINKING
INCLUDES+=-I $(CUDA_PATH)include/
LIB_PATH+=-L $(CUDA_LIB_PATH)


DYNAMIC_LIBS=-Wl,-rpath,"$(AF_PATH)lib"
####
####If you have problems linking with CUDA, try uncommenting the next line
####
DYNAMIC_LIBS+=-Wl,-rpath,"$(CUDA_LIB_PATH)"

LINKS = $(INCLUDES) $(DYNAMIC_LIBS) $(LIB_PATH) $(LIBS) 

objects = main.o

SHELL := /bin/bash

gpu : $(objects)
		g++ $(CFLAGS) -o ftle_gpu $(objects) $(LINKS) 
cpu : $(objects)
		g++ $(CFLAGS) -o ftle_cpu $(objects) $(LINKS) -lafcpu 

main.o: main.cpp  initializing.h initializing_inline.h \
		part_trans.h part_trans_inline.h testing.h testing_inline.h
		g++ $(CFLAGS) -c -o  main.o main.cpp $(LINKS)



remake : clean all

.PHONY : clean
clean :
	rm -f *.o && rm -f *.gch && rm -f && rm -f *.so* && rm ftle_*pu

.PHONY : run
run : clean omp
	./meme.out
