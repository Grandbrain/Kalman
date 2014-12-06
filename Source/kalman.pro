#-------------------------------------------------
#
# Project created by QtCreator 2014-12-06T18:27:05
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = kalman
TEMPLATE = app


SOURCES += main.cpp\
        mainwindow.cpp

HEADERS  += mainwindow.h \
    basic.h \
    extended.h \
    matrix.h \
    vector.h

FORMS    += mainwindow.ui
