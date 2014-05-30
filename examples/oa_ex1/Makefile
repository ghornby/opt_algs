PROJECT_ROOT_PATH = $(HOME)

COMMON_CLASSES_VERSION = cmnclass.1
COMMON_CLASSES_LIB_NAME = lib$(COMMON_CLASSES_VERSION).a
COMMON_CLASSES_INC_PATH = $(PROJECT_ROOT_PATH)/include/$(COMMON_CLASSES_VERSION)/

OPTALG_VERSION = optalg.1
OPTALG_LIB_NAME = lib$(OPTALG_VERSION).a
OPTALG_INC_PATH = $(PROJECT_ROOT_PATH)/include/$(OPTALG_VERSION)/

PROJECT_INCLUDES = -I$(PROJECT_ROOT_PATH)/include/ \
		-I$(COMMON_CLASSES_INC_PATH) \
		-I$(OPTALG_INC_PATH)


#OPTFLAGS = -Wall -O3 -march=i686 -fforce-addr -fomit-frame-pointer -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
#OPTFLAGS = -Wall -O2 -march=i686 -fforce-addr -fomit-frame-pointer -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
#OPTFLAGS = -Wall -O2 -fforce-addr -fomit-frame-pointer -funroll-loops -frerun-cse-after-loop -frerun-loop-opt -falign-functions=4
OPTFLAGS = -g -Wall


CFLAGS = $(OPTFLAGS) $(PROJECT_INCLUDES)


CC = g++



LINK =  \
	-l$(COMMON_CLASSES_VERSION) \
	-l$(OPTALG_VERSION) \
	-L$(PROJECT_ROOT_PATH)/lib -lcmnclass \
	-lm -lstdc++ -lgmp -lgmpxx \
	$(MYSQL_LINK) $(OPENGL_LINK) -pthread



TARGETS = \
	oa_ex1


all:	$(TARGETS)


clean:
	$(RM) *.o
	$(RM) $(TARGETS)


main.o:		main.cpp
		$(CC) $(CFLAGS) -c main.cpp

oa_ex1:		main.o \
		$(PROJECT_ROOT_PATH)/lib/lib$(OPTALG_VERSION).a
		$(CC) $(CFLAGS) main.o $(LINK) -o oa_ex1

