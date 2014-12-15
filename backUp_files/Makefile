ROOT = $(realpath ./)
EXT_PATH=$(realpath external)
USE_CPPITERTOOLS = 1
USE_CPPPROGUTILS = 1
USE_R = 1
USE_BOOST = 1
USE_ZI_LIB = 1
USE_BAMTOOLS = 1
USE_CURL = 1

include $(ROOT)/makefile-common.mk

UNAME_S := $(shell uname -s)
# from http://stackoverflow.com/a/8654800
OBJ_DIR = build
OBJ = $(addprefix $(OBJ_DIR)/, $(patsubst %.cpp, %.o, $(call rwildcard, src/, *.cpp)))
HEADERS = $(call rwildcard, src/, *.h) \
	$(call rwildcard, src/, *.hpp)

BIN = bin/sequenceTools
BIOALG = bin/bioalg
PROTO = bin/proto
EULER = bin/euler
QTTEST = bin/qtTest

COMMON = $(COMMON_OPT) -I/Users/nickhathaway/Documents/compilerProfiler/src
#COMMON = $(COMMON_DEBUG)
LD_FLAGS +=-L/Users/nickhathaway/Documents/compilerProfiler/lib -Wl,-rpath,/Users/nickhathaway/Documents/compilerProfiler/lib -lmacClangProfiler
############ main
.PHONY: all
all: $(OBJ_DIR) $(BIN)  
	scripts/fixDyLinking_mac.sh bin external

############ main
.PHONY: dev
dev: $(OBJ_DIR) $(BIN) $(PROTO) $(BIOALG) $(EULER) $(QTTEST) 
	scripts/fixDyLinking_mac.sh bin external

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)
	mkdir -p bin

$(OBJ_DIR)/%.o: %.cpp $(HEADERS)
	@mkdir -p $(OBJ_DIR)/$(shell dirname $<)
	$(CPP) -DNOT_HEADER_ONLY $(COMMON) -c $< -o $@

$(OBJ_DIR)/%.d: %.cpp
	$(SHELL) -ec '$(CPP) -M  $< | sed '\"s/$*.o/& $@/g'\" > $@'

$(BIN): $(OBJ) build/src/main.o
	$(CPP) -o $@ $^ $(LD_FLAGS)

$(PROTO): $(OBJ) build/src/main.o
	$(CPP) -o $@ $^ $(LD_FLAGS) 
$(BIOALG): $(OBJ) build/src/main.o
	$(CPP) -o $@ $^ $(LD_FLAGS)
	
$(EULER): $(OBJ) build/src/main.o
	$(CPP) -o $@ $^ $(LD_FLAGS)
	
$(QTTEST): $(OBJ) build/src/main.o
	$(CPP) -o $@ $^ $(LD_FLAGS)
############ clean
.PHONY: clean
clean:
	@rm -f $(BIN) $(PROTO) $(BIOALG) $(EULER) $(QTTEST)
	@rm -rf $(OBJ_DIR)

############ header
.PHONY: header
header:
	$(CPP) -o $(BIN) src/main.cpp $(COMMON) $(LD_FLAGS)

############ profileHeader
.PHONY: profileHeader
profileHeader:
	$(CPP) -o $(BIN)_profile src/main.cpp -pg $(COMMON) $(LD_FLAGS)
