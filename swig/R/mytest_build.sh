#!/bin/bash

# build for test
./clean
./docmatlab2R r_spams
./mybuild -nb
R CMD check pkg/spams
