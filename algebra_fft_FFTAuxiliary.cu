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
#include <gmp.h>
#include "cgbn/cgbn.h"

using namespace std;



#define CUDA_CALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl;  \
}

class fft_params_t {
    public:
    // parameters used by the CGBN context
    static const uint32_t TPB=0;                     // get TPB from blockDim.x  
    static const uint32_t MAX_ROTATION=4;            // good default value
    static const uint32_t SHM_LIMIT=0;               // no shared mem available
    static const bool     CONSTANT_TIME=false;       // constant time implementations aren't available yet
    
    // parameters used locally in the application
    static const uint32_t TPI=4;                   // threads per instance
    static const uint32_t BITS=512;                 // instance size
    static const uint32_t num_of_bytes=64;                 // instance size

};

typedef cgbn_mem_t<fft_params_t::BITS> Scalar;
typedef cgbn_context_t<fft_params_t::TPI>   context_t;
typedef cgbn_env_t<context_t, fft_params_t::BITS>   env_t;
typedef typename env_t::cgbn_t                 bn_t;
typedef typename env_t::cgbn_local_t          bn_local_t;



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


__global__ void cuda_fft_first_step( Scalar *input_field, Scalar omega, const size_t length, const size_t log_m) {

    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/fft_params_t::TPI;
    //printf("blockIdx=%d, blockDim.x=%d, threadIdx.x=%d, idx=%d\n",blockIdx.x, blockDim.x, threadIdx.x, idx);
    //const size_t block_length =( 1ul <<  LOG_NUM_THREADS) / fft_params_t::TPI;
    if(idx > length)
        return;
    context_t _context;
    env_t    _env(_context);

    //printf("CGBN with idx=%d, block_length=%d, startidx=%d\n", idx, block_length, startidx);
    bn_t  a, b; 
    /* swapping in place (from Storer's book) */
    size_t global_k = idx;
    size_t rk = bitreverse(global_k, log_m);
    if (global_k < rk  && rk < length)
    {
        cgbn_load(_env,a, &(input_field[global_k]));
        cgbn_load(_env,b, &(input_field[rk]));
        cgbn_store(_env, &(input_field[global_k]), b);
        cgbn_store(_env, &(input_field[rk]), a);
    }
    
   // __syncthreads();
    
}


__global__ void cuda_fft_second_step(Scalar *input_field, Scalar omega_binary, const size_t length, const size_t log_m, size_t s_index) {
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/fft_params_t::TPI;

    context_t _context;
    env_t    _env(_context);

    Scalar modulus_binary;
    //Fr modulus is : |811880050|3778125865|3092268470|2172737629|674490440|2042196113|1138881939|4026531841|
    uint32_t modulus_raw[16] = {4026531841,1138881939,2042196113,674490440, 
                                2172737629,3092268470,3778125865,811880050,
                                0, 0,0,0,
                                0,0,0,0};
    
    
//lianke: below is the modulus for FakeProofTest
// {1,6144,0,0, 
//                                 0,1048576,0,0,
//                                 0, 0,0,0,
//                                 0,0,0,0};
    memcpy(modulus_binary._limbs, modulus_raw, fft_params_t::num_of_bytes);
    bn_t modulus;
    cgbn_load(_env, modulus, &modulus_binary);
    bn_t omega;
    cgbn_load(_env, omega, &omega_binary);

    size_t m = 1 << (s_index - 1); // invariant: m = 2^{s-1}
    Scalar exponential_binary;
    exponential_binary._limbs[0] = (uint32_t)length/(2*m);
    bn_t exponential, w_m;
    cgbn_load(_env, exponential, &exponential_binary);
    cgbn_modular_power(_env, w_m, omega, exponential, modulus);
    // w_m is 2^s-th root of unity now

    size_t global_k = (idx / m) * m * 2 + idx %  m;
    if(global_k < length){
        //printf("global_k=%d\n", global_k);
        bn_t w, w_exp;
        size_t w_exp_int = idx % m;
        Scalar w_exp_binary;
        w_exp_binary._limbs[0] = w_exp_int;
        cgbn_load(_env, w_exp, &w_exp_binary);
        cgbn_modular_power(_env, w, w_m, w_exp, modulus);

        bn_t t;
        bn_t input_kjm, input_kj;
        cgbn_load(_env, input_kjm, &input_field[global_k +m]);
        cgbn_load(_env, input_kj, &input_field[global_k]);
        cgbn_mul(_env, t, w, input_kjm);
        cgbn_rem(_env, t, t, modulus);

        cgbn_add(_env, input_kj, input_kj, modulus);
        cgbn_sub(_env, input_kj, input_kj, t);
        //after subtraction, the result could be negative. so i need to add modulus to make it a positive number 
        cgbn_rem(_env, input_kj, input_kj, modulus);

        cgbn_store(_env, &input_field[global_k + m], input_kj);
        cgbn_load(_env, input_kj, &input_field[global_k]);
        cgbn_add(_env, input_kj, input_kj, t);
        cgbn_rem(_env, input_kj, input_kj, modulus);

        cgbn_store(_env, &input_field[global_k], input_kj);
        cgbn_mul(_env, w, w, w_m);
        cgbn_rem(_env, w, w, modulus);
    }
   // __syncthreads();

}




void best_fft (std::vector<Scalar> &a, const Scalar &omg)
{
	int cnt;
    cudaGetDeviceCount(&cnt);

    printf("CUDA Devices: %d, input_field size: %lu, input_field count: %lu\n", cnt, sizeof(Scalar), a.size());

    size_t threads_per_block = 128;
    size_t instance_per_block = (threads_per_block / fft_params_t::TPI);//TPI threads per instance, each block has threads.
    size_t blocks = (a.size() + instance_per_block - 1) / instance_per_block;

    printf("num of blocks %lu, threads per block %lu \n", blocks, threads_per_block);
    CUDA_CALL(cudaSetDevice(0));
    Scalar *in; 
    CUDA_CALL( cudaMalloc((void**)&in, sizeof(Scalar) * a.size()); )
    CUDA_CALL( cudaMemcpy(in, (void**)&a[0], sizeof(Scalar) * a.size(), cudaMemcpyHostToDevice); )

    const size_t length = a.size();
    const size_t log_m = log2(length); 
    //auto start = std::chrono::steady_clock::now();
    printf("launch block = %d thread = %d\n", blocks, threads_per_block);
    cuda_fft_first_step <<<blocks,threads_per_block>>>( in, omg, length, log_m);
    CUDA_CALL(cudaDeviceSynchronize());
    //cout << "finish first round" <<endl;
    size_t s = 1;
    for(; s <= log_m; s++){
        cuda_fft_second_step <<<blocks,threads_per_block>>>( in, omg, length, log_m, s);
        CUDA_CALL(cudaDeviceSynchronize());
        //cout <<"finish round " << s  <<endl;
    }

    // auto end = std::chrono::steady_clock::now();
    // std::chrono::duration<double> elapsed_seconds = end-start;
    // std::cout << "CUDA FFT elapsed time: " << elapsed_seconds.count() << "s\n";

    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }

    CUDA_CALL(cudaMemcpy((void**)&a[0], in, sizeof(Scalar) * a.size(), cudaMemcpyDeviceToHost); )
    CUDA_CALL(cudaDeviceSynchronize());
    CUDA_CALL(cudaFree(in));
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
        char* tmp = (char*)&inputArray[i]._limbs;

        memcpy(tmp, bytes, len);
    }


    Scalar omega;
    char* bytes = (char*)env->GetByteArrayElements(omegaArray, NULL);
    int len = env->GetArrayLength(omegaArray);
    char* tmp = (char*)&omega._limbs;
    memcpy(tmp , 
                bytes,
                len);
    
    best_fft(inputArray, omega);


    jbyteArray resultByteArray = env->NewByteArray((jsize)fft_params_t::num_of_bytes * input_len);
    for(int i=0; i < input_len;i++){
        // cout <<"cpp side output=";
        // inputArray[i].printBinaryDebug();
        env->SetByteArrayRegion(resultByteArray, i * fft_params_t::num_of_bytes , fft_params_t::num_of_bytes,   reinterpret_cast< jbyte*>(&inputArray[i]._limbs));
    }


    return resultByteArray;

}

