#-------------------------------------------------
#
# Project created by QtCreator 2017-07-10T21:09:21
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = Finance_Lab
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp \
    header_files/exotic_options.cpp \
    header_files/integration_functions.cpp \
    header_files/multivariate_integration.cpp \
    header_files/random_functions.cpp \
    header_files/simulation_functions.cpp

HEADERS  += mainwindow.h \
    header_files/exotic_options.hpp \
    header_files/integration_functions.hpp \
    header_files/multivariate_integration.hpp \
    header_files/random_functions.hpp \
    header_files/simulation_functions.hpp

FORMS    += mainwindow.ui

CONFIG += c++11

LIBS += -LC:/gsl/lib -lgsl -lgslcblas -lm
