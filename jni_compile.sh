JAVA_HOME='/usr/lib/jvm/java-8-openjdk-amd64'
NVCC_FLAGS :=-O3 -I include

sudo rm lib*.so

g++ -o3 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC  algebra_msm_FixedBaseMSM.cc -o libAlgebraMSMFixedBaseMSM.so 
g++ -o3 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC  algebra_msm_VariableBaseMSM.cc -o libAlgebraMSMVariableBaseMSM.so 
g++ -o3 -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC  algebra_fft_FFTAuxiliary.cc -o libAlgebraFFTAuxiliary.so 

sudo cp lib*.so /usr/lib/
javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/  src/main/java/algebra/msm/FixedBaseMSM.java 
javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/  src/main/java/algebra/msm/VariableBaseMSM.java 
javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/  src/main/java/algebra/msm/VariableBaseMSM.java 

cp src/main/java/algebra/msm/FixedBaseMSM.class target/classes/algebra/msm/
cp src/main/java/algebra/msm/VariableBaseMSM.class target/classes/algebra/msm/
cp src/main/java/algebra/fft/FFTAuxiliary.class target/classes/algebra/fft/

#compile the whole project and skip unit test.
#mvn install -DskipTests

#run a java class
#CLASSPATH=target/classes/:target/test-classes:dependency/spark-core_2.10-2.2.2.jar:dependency/scala-library-2.10.6.jar java algebra.msm.FixedBaseMSM 

#run a java test class 
# CLASSPATH=target/classes/:target/test-classes:dependency/junit-4.11.jar java junit.textui.TestRunner algebra.msm.SerialFixedBaseMSMTest