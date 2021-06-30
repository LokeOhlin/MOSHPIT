#/bin/bash

# Some general parameters
Emin="--Emin 8.9712e-12"
Emax="--Emax 3.204e-11"
Bins="--numBinsSubIon 1 --numBinsFullIon 1"

Eflags="${Emin} ${Emax} ${Bins}" 

workdir=$(pwd)
# Move to compilation directory, set correct compilation flags and then compile
cd ../src/bin
cp PDR_Makefile Makefile
#cp chem_Makefile Makefile
make clean
make
cd ${workdir}

# test V1
cd V1
pwd
/usr/bin/python3 ../radiationFields.py ${Eflags} --dxOrt 2.44e+14 --strength 10
cp ../../src/bin/moshpit .
./moshpit
cd ../

# test V2
cd V2
/usr/bin/python3 ../radiationFields.py ${Eflags} --dxOrt 2.44e+14 --strength 1e5
cp ../../src/bin/moshpit .
./moshpit
cd ../

# test V3
cd V3
/usr/bin/python3 ../radiationFields.py ${Eflags} --dxOrt 6.1e+8 --strength 10
cp ../../src/bin/moshpit .
./moshpit
cd ../

# test V4
cd V4
/usr/bin/python3 ../radiationFields.py ${Eflags} --dxOrt 2.44e+12 --strength 1e5
cp ../../src/bin/moshpit .
./moshpit
cd ../
