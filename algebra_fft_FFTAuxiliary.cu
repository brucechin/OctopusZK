#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <cmath>
#include <chrono>
#include "algebra_fft_FFTAuxiliary.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include "device_field.h"

using namespace std;

#define LOG_NUM_THREADS 10
#define NUM_THREADS (1 << LOG_NUM_THREADS)
#define LOG_CONSTRAINTS 20
#define CONSTRAINTS (1 << LOG_CONSTRAINTS)

#define CUDA_CALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl;  \
}

int reverseBits(int n, int range) {
    int ans = 0;
    for(int i = range - 1; i >= 0; i--){
        ans |= (n & 1) <<i;
        n>>=1;
    }
    return ans;
}


__device__ __forceinline__
size_t bitreverse(size_t n, const size_t l)
{
    return __brevll(n) >> (64ull - l); 
}

// __device__ uint32_t _mod [SIZE] = {0,0,0,
//                             0,0,0,0,
//                             0,0,1048576,
//                             0,0,0,6144,1};



__global__ void cuda_fft(Scalar *out, Scalar *field, const size_t length, const size_t log_m) {
    const int idx = blockIdx.x * blockDim.x + threadIdx.x;
    //printf("blockIdx=%d, blockDim.x=%d, threadIdx.x=%d, idx=%d\n",blockIdx.x, blockDim.x, threadIdx.x, idx);
    const size_t block_length = 1ul << (LOG_CONSTRAINTS - LOG_NUM_THREADS) ;
    const size_t startidx = idx * block_length;
    if(startidx > length)
        return;
    Scalar a[block_length];
    //Scalar* a = new Scalar[block_length];
    //TODO algorithm is non-deterministic because of padding
    Scalar omega_j = Scalar(_mod);
    omega_j = omega_j ^ idx; // pow
    Scalar omega_step = Scalar(_mod);
    omega_step = omega_step ^ (idx << (LOG_CONSTRAINTS - LOG_NUM_THREADS));
    
    Scalar elt = Scalar::one();
    //Do not remove log2f(n), otherwise register overflow
    size_t n = block_length, logn = log2f(n);
    assert (n == (1u << logn));
    for (size_t i = 0; i < 1ul<<(LOG_CONSTRAINTS - LOG_NUM_THREADS); ++i)
    {
        const size_t ri = bitreverse(i, logn);
        for (size_t s = 0; s < NUM_THREADS; ++s)
        {
            // invariant: elt is omega^(j*idx)
            size_t id = (i + (s<<(LOG_CONSTRAINTS - LOG_NUM_THREADS))) % (1u << log_m);
            Scalar tmp = field[id];
            tmp = tmp * elt;
            if (s != 0) tmp = tmp + a[ri];
            a[ri] = tmp;
            elt = elt * omega_step;
        }
        elt = elt * omega_j;
    }

    const Scalar omega_num_cpus = Scalar(_mod) ^ NUM_THREADS;
    size_t m = 1; // invariant: m = 2^{s-1}
    for (size_t s = 1; s <= logn; ++s)
    {
        // w_m is 2^s-th root of unity now
        const Scalar w_m = omega_num_cpus^(n/(2*m));
        for (size_t k = 0; k < n; k += 2*m)
        {
            Scalar w = Scalar::one();
            for (size_t j = 0; j < m; ++j)
            {
                const Scalar t = w;
                w = w * a[k+j+m];
                a[k+j+m] = a[k+j] - t;
                a[k+j] = a[k+j] + t;
                w = w * w_m;
            }
        }
        m = m << 1;
    }
    for (size_t j = 0; j < 1ul<<(log_m - LOG_NUM_THREADS); ++j)
    {
        if(((j << LOG_NUM_THREADS) + idx) < length)
            out[(j<<LOG_NUM_THREADS) + idx] = a[j];
    }
    
}




void best_fft (std::vector<Scalar> &a, const Scalar &omg)
{
	int cnt;
    cudaGetDeviceCount(&cnt);
    printf("CUDA Devices: %d, Field size: %lu, Field count: %lu\n", cnt, sizeof(Scalar), a.size());

    size_t blocks = NUM_THREADS / 256 + 1;
    size_t threads = NUM_THREADS > 256 ? 256 : NUM_THREADS;
    printf("NUM_THREADS %u, blocks %lu, threads %lu \n",NUM_THREADS, blocks, threads);

    Scalar *in;
    CUDA_CALL( cudaMalloc((void**)&in, sizeof(Scalar) * a.size()); )
    CUDA_CALL( cudaMemcpy(in, (void**)&a[0], sizeof(Scalar) * a.size(), cudaMemcpyHostToDevice); )

    Scalar *out;
    CUDA_CALL( cudaMalloc(&out, sizeof(Scalar) * a.size()); )
    const size_t length = a.size();
    const size_t log_m = log2(length); 
    //auto start = std::chrono::steady_clock::now();
    cuda_fft <<<blocks,threads>>>(out, in, length, log_m);
    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // std::cout << "CUDA FFT elapsed time: " << elapsed_seconds.count() << "s\n";

    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }

    CUDA_CALL( cudaMemcpy((void**)&a[0], out, sizeof(Scalar) * a.size(), cudaMemcpyDeviceToHost); )

    CUDA_CALL( cudaDeviceSynchronize();)
}

/*
 * Class:     algebra_fft_FFTAuxiliary
 * Method:    serialRadix2FFTNativeHelper
 * Signature: (Ljava/util/List;[B)[B
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_fft_FFTAuxiliary_serialRadix2FFTNativeHelper
  (JNIEnv * env, jclass obj, jobject inputs, jbyteArray omegaArray){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");
    jint input_len = env->CallIntMethod(inputs, java_util_ArrayList_size);


    vector<Scalar> inputArray = vector<Scalar>(input_len, Scalar());
    //TODO lianke update copy from java
    for(int i =0; i < input_len; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(inputs, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)&inputArray[i].im_rep;

        memcpy(tmp +Scalar::num_of_bytes - len, 
                                    bytes,
                                    len);
    }


    Scalar omega;
    char* bytes = (char*)env->GetByteArrayElements(omegaArray, NULL);
    int len = env->GetArrayLength(omegaArray);
    char* tmp = (char*)&omega.im_rep;
    memcpy(tmp +Scalar::num_of_bytes - len, 
                            bytes,
                            len);
    
    best_fft(inputArray, omega);


    jbyteArray resultByteArray = env->NewByteArray((jsize)Scalar::num_of_bytes * input_len);

    for(int i=0; i < input_len;i++){
        // cout <<"cpp side output=";
        // inputArray[i].printBinaryDebug();
        env->SetByteArrayRegion(resultByteArray, i * Scalar::num_of_bytes , Scalar::num_of_bytes,   reinterpret_cast<const jbyte*>(inputArray[i].im_rep));
    }


    return resultByteArray;

}

