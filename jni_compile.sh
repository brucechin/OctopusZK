JAVA_HOME='/usr/lib/jvm/java-8-openjdk-amd64'

sudo rm lib*.so
g++ -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC  algebra_msm_FixedBaseMSM.cc -o libAlgebraMSMFixedBaseMSM.so 
g++ -I $JAVA_HOME/include -I $JAVA_HOME/include/linux -shared -fPIC  algebra_msm_VariableBaseMSM.cc -o libAlgebraMSMVariableBaseMSM.so 
sudo cp lib*.so /usr/lib/
javac  -cp /home/lianke/dizk/dependency/scala-library-2.10.6.jar:/home/lianke/dizk/dependency/spark-core_2.10-2.2.2.jar:/home/lianke/dizk/target/classes/  src/main/java/algebra/msm/FixedBaseMSM.java 
javac  -cp /home/lianke/dizk/dependency/scala-library-2.10.6.jar:/home/lianke/dizk/dependency/spark-core_2.10-2.2.2.jar:/home/lianke/dizk/target/classes/  src/main/java/algebra/msm/VariableBaseMSM.java 
cp src/main/java/algebra/msm/FixedBaseMSM.class target/classes/algebra/msm/
cp src/main/java/algebra/msm/VariableBaseMSM.class target/classes/algebra/msm/


#run a java class
#CLASSPATH=target/classes/:target/test-classes:dependency/spark-core_2.10-2.2.2.jar:dependency/scala-library-2.10.6.jar java algebra.msm.FixedBaseMSM 

#run a java test class 
# CLASSPATH=target/classes/:target/test-classes:dependency/junit-4.11.jar java junit.textui.TestRunner algebra.msm.SerialFixedBaseMSMTest