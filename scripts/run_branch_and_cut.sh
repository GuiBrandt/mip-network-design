#!/usr/bin/env bash

mkdir -p logs

for size in 10 20 50
do
for n in 1 2 3 4 5
do
for instance in "instances/mo420_network_design_${n}_${size}_"*
do
    echo $size $n $instance
    time build/mo420-network-design-2 $instance > logs/branch-and-cut/$(basename $instance).log
done
done
done
