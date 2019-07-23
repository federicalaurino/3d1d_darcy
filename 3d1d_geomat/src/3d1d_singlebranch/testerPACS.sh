#!/bin/bash

inputfile=input.param
run="./M3D1D $inputfile"

echo "========================================"
echo "1 simulation SAMG direct solver"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[20,20,20]\';/" $inputfile
resu_folder=20sub_direct
printfile=20sub_direct.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk

echo "========================================"
echo "2 simulation SAMG direct solver"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[40,40,40]\';/" $inputfile
resu_folder=40sub_direct
printfile=40sub_direct.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk


echo "========================================"
echo "3 simulation SAMG direct solver"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[60,60,60]\';/" $inputfile
resu_folder=60sub_direct
printfile=60sub_direct.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk


echo "========================================"
echo "4 simulation SAMG direct solver"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[80,80,80]\';/" $inputfile
resu_folder=80sub_direct
printfile=80sub_direct.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk

