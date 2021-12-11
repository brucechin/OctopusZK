# 2021.12.8

Implement the JNI code for `FixedBaseMSM.batchMSM` and `FixedBaseMSM.doubleBatchMSM` on cpp side. Reduce time consumption over naive Java DIZK version by 20% on microbenchmarks. Next step to migrate these two APIs to CUDA code.


# 2021.12.10

Implement the correct cpp logic for FixedBaseMSM and return to java side for correctness check for fakeG1 and fakeG2. Do some profiling and find that the data movement across java and cpp brings little overhead compared to the heavy computation on large finite field and curves. next step is to implement BN254 curves. 

Start a single node spark service and successfully submit the jobs to it. Profiling will not be a problem.

Look into how to design a GPU-friendly BigInt library.