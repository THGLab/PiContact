#! /bin/bash

SUBMIT=$1

for TZVPD in *_TZVPD.out
do

NAME=${TZVPD%.out}

echo $NAME >>name.dat
#echo $NAME >>elec.dat
grep  "E_cls_elec" ${NAME}.out | awk '{print  $6}' >>elec.dat
#echo $NAME >>pauli.dat
grep  "E_mod_pauli" ${NAME}.out | awk '{print  $6}' >>pauli.dat
#echo $NAME >>forzen.dat
grep  "FROZEN" ${NAME}.out | awk '{print  $2}' >>frozen.dat
#echo $NAME >>dispersion.dat
grep  "DISPERSION" ${NAME}.out | awk '{print  $2}' >>dispersion.dat
#echo $NAME >>polarization.dat
grep  "POLARIZATION" ${NAME}.out | awk '{print  $2}' >>polarization.dat
#echo $NAME >>CT.dat
grep  "CHARGE TRANSFER" ${NAME}.out | awk '{print  $3}' >>CT.dat
#echo $NAME >>Tot.dat
grep  "TOTAL          " ${NAME}.out | awk '{print  $2}' >>Tot.dat

done

