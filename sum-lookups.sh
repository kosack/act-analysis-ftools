#!/bin/sh



# by telescope:
for val in MC_ENERGY HIL_TEL_WIDTH HIL_TEL_LENGTH; do
    for tel in 001 002 003 004; do
	for what in sum sum2 count; do
 	    echo SUMMING CT$tel $val $what ...
 	    python ${TOOLSDIR}/sum-maps.py -v -o CT${tel}-${val}-lookup-${what}.fits *-CT${tel}*${val}-lookup-${what}.fits
 	done
     done
 done


# by type
# for ttype in TYPE01_00; do 
#     for val in MC_ENERGY HIL_TEL_WIDTH HIL_TEL_LENGTH ; do
# 	for what in sum sum2 count; do
# 	    echo SUMMING $ttype $val $what ...
# 	    ${TOOLSDIR}/sum_maps.pl -o ${ttype}-${val}-lookup-${what}.fits \
# 		*-CT*-${ttype}-${val}-lookup-${what}.fits
# 	done
#     done
# done
