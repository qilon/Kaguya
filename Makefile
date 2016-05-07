# enable pretty build (comment to see full commands)
Q ?= @

# Set DEBUG to 1 for debugging
DEBUG := 0

ifeq ($(DEBUG), 1)
	BUILD_DIR := build_debug
else
	BUILD_DIR := build
endif

BIN_BUILD_DIR := $(BUILD_DIR)/bin

# All the directories containing code
SRC_DIRS := Kaguya

ALL_BUILD_DIRS := $(addprefix $(BUILD_DIR)/, $(SRC_DIRS))
ALL_BUILD_DIRS += $(LIB_BUILD_DIR)

EIGEN_INCLUDE := -I/usr/include/eigen3

FLAGS_INCLUDE := $(EIGEN_INCLUDE) -I./Kaguya 

# Library dependencies
GL_LIB := -lGL -lGLU -lglut -lX11 -lGLEW

BOOST_LIB := -lboost_filesystem -lboost_system -lboost_thread

OPENCV_LIB := -lopencv_core -lopencv_highgui -lopencv_imgproc -lopencv_imgcodecs

OPENMESH_LIB := -lOpenMeshCore -lOpenMeshTools

CERES_LIB := -lceres -lglog -ltbb -ltbbmalloc -lcholmod -lccolamd \
	-lcamd -lcolamd -lamd -lsuitesparseconfig -llapack -lf77blas -latlas

CLUSTER_LIB_DIR := /usr/local/lib

LIBRARY_DIRS := $(CLUSTER_LIB_DIR)

LDFLAGS := $(BOOST_LIB) $(OPENCV_LIB) $(GL_LIB) $(OPENMESH_LIB) $(CERES_LIB)
LDFLAGS += $(foreach library_dir, $(LIBRARY_DIRS), -L$(library_dir))
LDFLAGS += $(foreach library_dir, $(LIBRARY_DIRS), -Wl,-rpath,$(library_dir))

# Setting compiler and building flags
CXX := g++
CXXFLAGS += -std=c++11 -fopenmp $(FLAGS_INCLUDE)

# Debugging
ifeq ($(DEBUG), 1)
	CXXFLAGS += -DDEBUG -g -O0
else
	CXXFLAGS += -DNDEBUG -O2 -ffast-math -Wno-unused-result
endif

# Automatic dependency generation
CXXFLAGS += -MMD -MP

# Get all source files
CONSOLE_APP_SRCS := ./Kaguya/main.cpp ./Kaguya/BlackBox.cpp

CONSOLE_BIN := $(BIN_BUILD_DIR)/Kaguya

.PHONY: all

all: $(BIN_BUILD_DIR)
	$(CXX) $(CONSOLE_APP_SRCS) -o $(CONSOLE_BIN) $(CXXFLAGS) $(LINKFLAGS) $(LDFLAGS)

$(BIN_BUILD_DIR):
	@ mkdir -p $@
