#!/bin/zsh
ifort fit.f90
if [ $? -eq 0 ]; then
    ./a.out
else
    echo "Compile error"
    exit 1
fi

if [ $? -eq 0 ]; then
    python plot.py
    python ato-shori.py
else
    echo "Python error"
    exit 1
fi
