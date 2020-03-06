#!bin/bash
#install RNAmotif and RNAfold
set -e

cd dependencies/
tar -xvzf rnamotif-3.0.7.tar.gz
cd rnamotif-3.0.7/
make
cp src/rnamotif src/rmprune ../../ 
cd ../
#export PATH=$PATH:$(pwd)/src/rnamotif
echo "RNAmotif done!"
echo "RNAfold installation begins"
tar -xvzf ViennaRNA-2.1.8.tar.gz
cd ViennaRNA-2.1.8/
./configure --prefix=$pwd 
make
cp Progs/RNAfold ../../
#export PATH=$PATH:$(pwd)/Progs
echo "RNAfold installation finished!"
