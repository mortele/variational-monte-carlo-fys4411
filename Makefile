CXX:=g++

CXX_FLAGS:=-O3 -Wall -Wextra

HEADERS:=$(wildcard include/*.h)

APP_NAME = ${app}
SOURCES:=$(wildcard src/*.cpp)

OBJECTS:=$(SOURCES:.cpp=.o)

$(APP_NAME) : $(OBJECTS)
	@mkdir -p $(dir bin/$@)
	$(CXX) $(APP_NAME).cpp $^ -o bin/$@.out $(CXX_FLAGS)

clean:
	rm -f src/*.o
	rm -f bin/$(APP_NAME).out