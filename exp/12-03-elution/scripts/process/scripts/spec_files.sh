#!/bin/bash
matches=( 'G166' 'CB660' 'mES' 'Celegans' 'DSL2' 'Dembryo' 'Urchin' 'Ax4' )
species=( 'Hs' 'Hs' 'Mm' 'Ce' 'Dm' 'Dm' 'Sp' 'Dd' )
for spec in ${species[@]}
do
	mkdir $spec
done
for (( i = 0 ; i < ${#matches[@]} ; i++ ))
do
	mv *${matches[$i]}* ${species[$i]}
	# yadda yadda
done
