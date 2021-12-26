# 2021.12.8

Implement the JNI code for `FixedBaseMSM.batchMSM` and `FixedBaseMSM.doubleBatchMSM` on cpp side. Reduce time consumption over naive Java DIZK version by 20% on microbenchmarks. Next step to migrate these two APIs to CUDA code.


# 2021.12.10

Implement the correct cpp logic for FixedBaseMSM and return to java side for correctness check for fakeG1 and fakeG2. Do some profiling and find that the data movement across java and cpp brings little overhead compared to the heavy computation on large finite field and curves. next step is to implement BN254 curves. 

Start a single node spark service and successfully submit the jobs to it. Profiling will not be a problem.

Look into how to design a GPU-friendly BigInt library.

# 2021.12.11

previously, i convert java.BigInteger into byteArray and pass it to cpp side. cpp side bigint use char[] to store it. everything goes well.
However, when i try to store bigint with int[] at cpp side, due to the little/big endian problem, they are different.
Therefore, i need to define some toIntArray from java, and pass the intarray representation of java.BigInteger to cpp side, and use int[] to store it.

# 2021.12.13

kind of stuck in building lightweighted cpp fix length cpp bigint library. implement a simple modulus function. implement and verify the correctness of JNI cpp side `VariableBaseMSM.serialMSM` and `VariableBaseMSM.doubleMSM`. next step is to implement the FFT on JNI cpp code side. Also notice that the BN254a curve is not implemented yet.

# 2021.12.16

bug fixed:
(1). support a correct multiplication, mod, and pow with mod.
(2). fix a bug when (smaller - bigger) mod modulus.
(3). fix a bug when shift offset is 32: for example: 11 >> 32 = 11, instead of 0. however, 11 >> 31 = 0;

# 2021.12.19 

Because the GPUSnark repo is not correct with the DIZK computation results, and my naive lightweight cpp bigint library has some corner bugs. I decided to migrate my code to the second version of CGBN developed by NVIDIA. I spend some time figuring out how to divide the tasks to threads in GPU and finished the FFT in CUDA. next step is to implement for FixedbaseMSM and VariablebaseMSM.


# 2021.12.22 

meet some difficult issues when implementing MSM using the NVIDIA CGBN library. each cgbn_t is operated by TPI threads, so when i tried to write back the computation results back to CPU memory, i can only write a partial component with each thread. need to carefully operate on global variables and local variables in order to obtain the correct CGBN computation results and memcpy them back to CPU.


# 2021.12.25 
found out that the distributed zk can not execute. because of the following error
### Caused by: java.lang.IllegalArgumentException: Invalid lambda deserialization at reductions.r1cs_to_qap.R1CStoQAPRDD.$deserializeLambda$(R1CStoQAPRDD.java:1)'
 
i think it is due to some Serializable overriding bugs by myself?


found a stupid microbenchmark error in SerialFixedbaseMSMG1, they did not use the correct group factory
