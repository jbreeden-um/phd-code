# Computation of Viability Domains Under Multiple Control Barrier Functions with Input Constraints

### Overview

The main files in this folder are:
 * [algorithm.c](algorithm.c) contains the viability domain search algorithm. This file assumes that two CBFs are given and attempts to compute a third (or fourth depending on the case).
 * [expert.c](expert.c) contains tools to check a solution developed either by [algorithm.c](algorithm.c) or by hand.
 * [setup.c](setup.c) contains tools to verify that a single CBF is valid. This amounts to checking that the choice of `mu` is valid.
 
See the [Makefile](Makefile) for how to build all three programs.

### Viability Search Algorithm

The viability search algorithm in [algorithm.c](algorithm.c) has three modes. One can switch between the modes by changing the `#define` statements in [verify.h](verify.h).
 * `_RUN_CASE_2_` is the main version and the one discussed in the paper.
 * `_RUN_CASE_3_` is a smarter and faster version of `_RUN_CASE_2_` that forces all the unviable points found to belong to a single cluster, so that only one new CBF is introduced instead of two new CBFs.
 * `_RUN_CASE_1_` searches a fewer number of points than the other two cases, and thus runs far faster, but requires additional conservatism on the control set for its results to be accurate. This case is primarily included for experimental and rapid-prototyping purposes. 