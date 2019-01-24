BDIR := build
ODIR := $(BDIR)/obj
HEADERDIRS := . $(shell find * -path $(BDIR) -prune -o -type d -print)

CXX := g++
CXXFLAGS := -stdlib=c++11 -g -Wall -Wextra -I$(HEADERDIRS)


$(info $(HEADERDIRS))
