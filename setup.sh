sudo apt-get update 
sudo apt-get install -y maven openjdk-8-jdk g++ lzip scala git m4
sudo apt install -y nvidia-driver-470 nvidia-cuda-toolkit

#install GMP
wget https://gmplib.org/download/gmp/gmp-6.2.1.tar.lz
lzip -d  gmp-6.2.1.tar.lz
tar  -xvf gmp-6.2.1.tar
cd gmp-6.2.1/
export GMP_HOME=/home/lianke/gmp-6.2.1
./configure --prefix=$GMP_HOME
make -j16
make install


wget https://downloads.apache.org/spark/spark-3.2.0/spark-3.2.0-bin-hadoop3.2.tgz
tar xvf spark-*
sudo mv spark-3.2.0-bin-hadoop3.2 /opt/spark
echo "export SPARK_HOME=/opt/spark" >> ~/.profile
echo "export PATH=$PATH:$SPARK_HOME/bin:$SPARK_HOME/sbin" >> ~/.profile
echo "export PYSPARK_PYTHON=/usr/bin/python3" >> ~/.profile
source .profile