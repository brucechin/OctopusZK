sudo apt-get update 
sudo apt-get install maven openjdk-8-jdk g++ lzip 
sudo apt install nvidia-driver-470 nvidia-cuda-toolkit

#install GMP
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
lzip -d  gmp-6.2.1.tar.lz
tar  -xvf gmp-6.2.1.tar
cd gmp-6.2.1/
export GMP_HOME=/home/lianke/gmp-6.2.1
./configure --prefix=$GMP_HOME
make -j16
make install