CXX:=g++

CXX_FLAGS:=-Wall -Wextra -O3

HEADERS:=$(wildcard include/*.h)

APP_NAME = ${app}
SOURCES:=$(wildcard src/*.cpp)

OBJECTS:=$(SOURCES:.cpp=.o)

$(APP_NAME) : $(OBJECTS)
	@mkdir -p $(dir bin/$@)
	$(CXX) -g $(APP_NAME).cpp $^ -o bin/$@.out $(CXX_FLAGS)

clean:
	rm -f src/*.o
	rm -f bin/$(APP_NAME).out