#!/bin/bash

## build python pkg

# clean
./clean

# Interface building for test :
./mkpy spams

# Tests:
python test_spams.py test_prox
