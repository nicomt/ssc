TARGET_EXEC ?= sixparsolve.js

CC=emcc
CXX=em++
CPPFLAGS=-O3
LDFLAGS=-s EXPORTED_FUNCTIONS=_sixparsolve,_free,_malloc -s EXPORTED_RUNTIME_METHODS=cwrap -s ALLOW_MEMORY_GROWTH=1 -s MODULARIZE=1 -s FILESYSTEM=0 -s STRICT=1 -s ENVIRONMENT=node -s WASM=1 -s WASM_ASYNC_COMPILATION=0 -s LLD_REPORT_UNDEFINED --no-entry -s EXPORT_NAME="SixParSolve"
BUILD_DIR ?= ./build
SRC_DIRS ?= ./src

SRCS := $(shell find $(SRC_DIRS) -name *.cpp)
OBJS := $(SRCS:%=$(BUILD_DIR)/%.o)
DEPS := $(OBJS:.o=.d)

$(BUILD_DIR)/$(TARGET_EXEC): $(OBJS)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS)

# c++ source
$(BUILD_DIR)/%.cpp.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c $< -o $@

.PHONY: clean

clean:
	$(RM) -r $(BUILD_DIR)

-include $(DEPS)

MKDIR_P ?= mkdir -p
