#!/bin/bash

## build python pkg

# clean
./clean

# Interface building for test :
./mkpy -py3 spams

# Tests:
python3 test_spams.py
