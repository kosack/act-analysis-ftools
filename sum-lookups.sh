#!/bin/sh

TOOLSDIR=~/Source/Working/PyFITSTools

for tel in 001 002 003 004; do
    for val in HIL_TEL_WIDTH HIL_TEL_LENGTH MC_ENERGY; do
	for what in sum sum2 count; do
	    echo SUMMING CT$tel $val $what ...
	    ${TOOLSDIR}/sum_maps.pl -o CT${tel}-${val}-lookup-${what}.fits *-CT${tel}*${val}-lookup-${what}.fits
	done
    done
done