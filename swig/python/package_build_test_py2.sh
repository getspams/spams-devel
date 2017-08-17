#!/bin/bash

## build python pkg

# clean
./clean

# Interface building for test :
./mkpy spams

# Tests:
python2 test_spams.py
