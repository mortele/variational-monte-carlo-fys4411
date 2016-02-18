TEMPLATE = app
CONFIG  += console c++11
CONFIG  -= app_bundle
CONFIG  -= qt
@LIBS   += -L/usr/lib@
LIBS    += -lblas -llapack
INCLUDEPATH += /usr/local/Cellar/armadillo/6.200.4/include
QMAKE_CXXFLAGS += -O3
SOURCES += main.cpp \
    system.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/harmonicoscillator.cpp \
    particle.cpp \
    WaveFunctions/wavefunction.cpp \
    InitialStates/initialstate.cpp \
    InitialStates/randomuniform.cpp \
    Math/random.cpp \
    sampler.cpp \
    WaveFunctions/simplegaussian.cpp \
    WaveFunctions/gaussian4.cpp \
    WaveFunctions/multiparticleho.cpp \
    WaveFunctions/multiparticlehointeracting.cpp \
    Hamiltonians/harmonicoscillatorinteracting.cpp

HEADERS += \
    system.h \
    Hamiltonians/hamiltonian.h \
    Hamiltonians/harmonicoscillator.h \
    particle.h \
    WaveFunctions/wavefunction.h \
    InitialStates/initialstate.h \
    InitialStates/randomuniform.h \
    Math/random.h \
    sampler.h \
    WaveFunctions/simplegaussian.h \
    WaveFunctions/gaussian4.h \
    WaveFunctions/multiparticleho.h \
    examples.h \
    WaveFunctions/multiparticlehointeracting.h \
    Hamiltonians/harmonicoscillatorinteracting.h

