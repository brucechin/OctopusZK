# 2021.12.8

Implement the JNI code for `FixedBaseMSM.batchMSM` and `FixedBaseMSM.doubleBatchMSM` on cpp side. Reduce time consumption over naive Java DIZK version by 20% on microbenchmarks. Next step to migrate these two APIs to CUDA code.