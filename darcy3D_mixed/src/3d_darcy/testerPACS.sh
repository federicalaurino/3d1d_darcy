#!/bin/bash

inputfile=input.param
run="./M3D1D $inputfile"

echo "========================================"
echo "1 simulation gmresNOprec"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[10,10,10]\';/" $inputfile
resu_folder=10sub_gmresNOprec
printfile=10sub_gmresNOprec.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk

echo "========================================"
echo "2 simulation gmresNOprec"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[15,15,15]\';/" $inputfile
resu_folder=15sub_gmresNOprec
printfile=15sub_gmresNOprec.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk


echo "========================================"
echo "3 simulation gmresNOprec"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[20,20,20]\';/" $inputfile
resu_folder=20sub_gmresNOprec
printfile=20sub_gmresNOprec.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk


echo "========================================"
echo "4 simulation gmresNOprec"
echo "========================================"
sed -r -i "s/^NSUBDIV_T.*/NSUBDIV_T = \'[25,25,25]\';/" $inputfile
resu_folder=25sub_gmresNOprec
printfile=25sub_gmresNOprec.log
make
$run > $printfile
rm -rf $resu_folder
mv vtk $resu_folder
mv $printfile $resu_folder
mkdir vtk

