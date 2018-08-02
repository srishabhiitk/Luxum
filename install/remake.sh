cd ..
rm -rf build
rm -f ../luxum
mkdir build
cd build
cmake .. -DCMAKE_C_COMPILER=mpicc -DCMAKE_CXX_COMPILER=mpic++
make
