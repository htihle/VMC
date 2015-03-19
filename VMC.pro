TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    metropolis.cpp \
    expectationvalues.cpp \
    wavefunction.cpp \
    orbital.cpp \
    ho0.cpp \
    ho1.cpp \
    ho2.cpp \
    onedimensionalslater.cpp \
    slater.cpp \
    h1s.cpp \
    h2s.cpp \
    hamiltonian.cpp \
    atom.cpp \
    h2p.cpp \
    heliumwavefunction.cpp \
    wfTest.cpp

HEADERS += \
    metropolis.h \
    expectationvalues.h \
    wavefunction.h \
    orbital.h \
    ho0.h \
    ho1.h \
    ho2.h \
    onedimensionalslater.h \
    slater.h \
    h1s.h \
    h2s.h \
    hamiltonian.h \
    atom.h \
    h2p.h \
    heliumwavefunction.h \
    wfTest.h

    LIBS+= -larmadillo -lblas -llapack
