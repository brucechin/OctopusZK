JAVA_HOME='/usr/lib/jvm/java-8-openjdk-amd64'
GMP_HOME='~/OctopusZK/gmp-6.2.1'
CUDA_HOME='/usr/local/cuda-11.1/bin'
GPU_ARCH_CODE=sm_60

sudo rm lib*.so

nvcc -O3  -Xcompiler -fPIC -shared -I $GMP_HOME/include -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -I ./include -arch=$GPU_ARCH_CODE algebra_msm_FixedBaseMSM.cu -o libAlgebraMSMFixedBaseMSM.so 
nvcc  -O3 -Xcompiler  -fPIC -shared -I $GMP_HOME/include -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -I ./include -arch=$GPU_ARCH_CODE algebra_msm_VariableBaseMSM.cu -o libAlgebraMSMVariableBaseMSM.so 
nvcc   -Xcompiler  -fPIC -shared -I $GMP_HOME/include -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -I ./include -arch=$GPU_ARCH_CODE algebra_fft_FFTAuxiliary.cu -o libAlgebraFFTAuxiliary.so

sudo cp lib*.so /usr/lib/

# javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/:./dependency/*  src/main/java/algebra/msm/FixedBaseMSM.java 
# javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/:./dependency/*  src/main/java/algebra/msm/VariableBaseMSM.java 
# javac  -cp ./dependency/org-apache-commons-codec.jar:./dependency/scala-library-2.10.6.jar:./dependency/spark-core_2.10-2.2.2.jar:./target/classes/:./dependency/*  src/main/java/algebra/fft/FFTAuxiliary.java 

# cp src/main/java/algebra/msm/FixedBaseMSM.class target/classes/algebra/msm/
# cp src/main/java/algebra/msm/VariableBaseMSM.class target/classes/algebra/msm/
# cp src/main/java/algebra/fft/FFTAuxiliary.class target/classes/algebra/fft/

#compile the whole project and skip unit test.
#mvn install -DskipTests

#run a java class
#CLASSPATH=target/classes/:target/test-classes:dependency/spark-core_2.10-2.2.2.jar:dependency/scala-library-2.10.6.jar java algebra.msm.FixedBaseMSM 
