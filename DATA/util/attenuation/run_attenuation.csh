#!/bin/csh

echo "computing constants for a range of periods"

rm -f q* Q* modulus*

xcompute < input_Q_DK_01_1000 > result01
xcompute < input_Q_DK_02_1000 > result02
xcompute < input_Q_DK_03_1000 > result03
xcompute < input_Q_DK_04_1000 > result04
xcompute < input_Q_DK_05_1000 > result05
xcompute < input_Q_DK_08_1000 > result08
xcompute < input_Q_DK_10_1000 > result10
xcompute < input_Q_DK_20_1000 > result20

rm -f q* Q* modulus*

