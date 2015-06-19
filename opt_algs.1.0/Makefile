OPTALG_NAME = optalg.1
OPTALG_VERSION = optalg.1.0
OPTALG_LIB_NAME = lib$(OPTALG_NAME).a
OPTALG_LIB_VERSION = lib$(OPTALG_VERSION).a


COMMON_CLASSES_VERSION = cmnclass.1
COMMON_CLASSES_LIB_NAME = lib$(COMMON_CLASSES_VERSION).a
COMMON_CLASSES_INC_PATH = $(HOME)/include/$(COMMON_CLASSES_VERSION)/


PROJECT_INCLUDES = -I$(HOME)/include/ \
		-I$(COMMON_CLASSES_INC_PATH)


#OPTFLAGS = -Wall -O3 -march=i686 -fforce-addr -fomit-frame-pointer -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
#OPTFLAGS = -Wall -O2
OPTFLAGS = -g -Wall
CFLAGS = -std=c++11 $(OPTFLAGS)
INCL =	$(PROJECT_INCLUDES)

UNAME := $(shell uname -s)
ifeq ($(UNAME), Darwin)
INCL := $(INCL) -I/opt/local/include
endif

LINK = -L$$(HOME)/lib -l$(COMMON_CLASSES_VERSION) -lm -lstdc++
ifeq ($(UNAME), Darwin)
LINK := $(LINK) -L/opt/local/lib
endif



CC = g++

all:	$(OPTALG_LIB_VERSION)


OPTALG_OBJS =	individual.o ind_history.o individ_real.o opt_alg.o \
		oa_random.o oa_hillclimb.o oa_function.o oa_frosen.o


TESTS = 

TESTDIR = tests/

$(OPTALG_LIB_VERSION): $(OPTALG_OBJS)
	ar rcv $(OPTALG_LIB_VERSION) $(OPTALG_OBJS)


opt_alg.o:	opt_alg.cpp opt_alg.h
		$(CC) $(CFLAGS) $(INCL) -c opt_alg.cpp


oa_hillclimb.o:	oa_hillclimb.cpp oa_hillclimb.h opt_alg.h
		$(CC) $(CFLAGS) $(INCL) -c oa_hillclimb.cpp


oa_random.o:	oa_random.cpp oa_random.h opt_alg.h
		$(CC) $(CFLAGS) $(INCL) -c oa_random.cpp



oa_function.o:	oa_function.cpp oa_function.h
		$(CC) $(CFLAGS) $(INCL) -c oa_function.cpp


oa_frosen.o:	oa_frosen.cpp oa_frosen.h oa_function.h
		$(CC) $(CFLAGS) $(INCL) -c oa_frosen.cpp







ind_history.o:	ind_history.cpp ind_history.h
		$(CC) $(CFLAGS) $(INCL) -c ind_history.cpp


individual.o:	individual.cpp individual.h ind_history.h
		$(CC) $(CFLAGS) $(INCL) -c individual.cpp


individ_real.o:	individ_real.cpp individ_real.h individual.h
		$(CC) $(CFLAGS) $(INCL) -c individ_real.cpp

#####################################################

GLOBAL_PATH = /usr/local

installglobal:
		cp $(OPTALG_LIB_VERSION) $(GLOBAL_PATH)/lib/$(OPTALG_LIB_VERSION)
		rm -f $(GLOBAL_PATH)/lib/$(OPTALG_LIB_NAME)
		ln -s $(OPTALG_LIB_VERSION) $(GLOBAL_PATH)/lib/$(OPTALG_LIB_NAME)
		rm -rf $(GLOBAL_PATH)/include/$(OPTALG_VERSION)
		-mkdir $(GLOBAL_PATH)/include/$(OPTALG_VERSION)
		cp *.h $(GLOBAL_PATH)/include/$(OPTALG_VERSION)
		rm -f $(GLOBAL_PATH)/include/$(OPTALG_NAME)
		ln -s $(OPTALG_VERSION) $(GLOBAL_PATH)/include/$(OPTALG_NAME)


installusr:
		cp $(OPTALG_LIB_VERSION) $(HOME)/lib/$(OPTALG_LIB_VERSION)
		rm -f $(HOME)/lib/$(OPTALG_LIB_NAME)
		ln -s $(OPTALG_LIB_VERSION) $(HOME)/lib/$(OPTALG_LIB_NAME)
		rm -rf $(HOME)/include/$(OPTALG_VERSION)
		-mkdir $(HOME)/include/$(OPTALG_VERSION)
		cp *.h $(HOME)/include/$(OPTALG_VERSION)
		rm -f $(HOME)/include/$(OPTALG_NAME)
		ln -s $(OPTALG_VERSION) $(HOME)/include/$(OPTALG_NAME)


.PHONY : clean

clean:
	$(RM) *.o
	$(RM) $(OPTALG_LIB_VERSION)


check: $(TESTS)
	cd $(TESTDIR) && ./$(TESTS)

