#!/usr/bin/env bash

N=8
T=8
N5=1
D=4

BETA_START=0.5
BETA_END=4.0
BETA_INC=0.5
EPS5=1.0

N_SWEEPS=20
N_DATA=20

./su2_x_test ${N} ${T} ${N5} ${D} ${BETA_START} ${BETA_END} ${BETA_INC} ${EPS5} ${N_SWEEPS} ${N_DATA}
