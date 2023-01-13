#!/bin/sh

for i in 1 2 3 4 5 6
do
	cd Melissa_Model$i
    ./3poppi.sh
	mv *3pop.pi.windowed.pi PI/3pop
    rm *3pop.pi.log
    cd ..
done