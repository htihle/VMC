TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    metropolis.cpp \
    expectationvalues.cpp \
    Orbitals/h1s.cpp \
    Orbitals/h2p.cpp \
    Orbitals/h2s.cpp \
    Orbitals/ho0.cpp \
    Orbitals/ho1.cpp \
    Orbitals/ho2.cpp \
    Orbitals/orbital.cpp \
    Hamiltonians/hamiltonian.cpp \
    Hamiltonians/atom.cpp \
    Wavefunctions/heliumwavefunction.cpp \
    Wavefunctions/onedimensionalslater.cpp \
    Wavefunctions/wavefunction.cpp \
    Wavefunctions/slater.cpp \
    Wavefunctions/wfTest.cpp \
    Hamiltonians/harmonicoscillator.cpp

HEADERS += \
    metropolis.h \
    expectationvalues.h \
    Orbitals/h1s.h \
    Orbitals/h2p.h \
    Orbitals/h2s.h \
    Orbitals/ho0.h \
    Orbitals/ho1.h \
    Orbitals/ho2.h \
    Orbitals/orbital.h \
    Hamiltonians/atom.h \
    Hamiltonians/hamiltonian.h \
    Wavefunctions/wfTest.h \
    Wavefunctions/wavefunction.h \
    Wavefunctions/slater.h \
    Wavefunctions/onedimensionalslater.h \
    Wavefunctions/heliumwavefunction.h \
    Hamiltonians/harmonicoscillator.h

    LIBS+= -larmadillo -lblas -llapack
