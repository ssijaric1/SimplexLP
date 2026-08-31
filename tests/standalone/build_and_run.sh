#!/usr/bin/env bash
# Builds the LP core against the shim headers and runs the test suite.
# No natID SDK required - pure standard C++20.
set -e
cd "$(dirname "$0")"
g++ -std=c++20 -O2 -Wall -Wextra -Ishims -o lp_tests test_main.cpp
./lp_tests
