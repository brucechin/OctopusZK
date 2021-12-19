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
//#include "device_field.h"
#include <gmp.h>
#include "cgbn/cgbn.h"

using namespace std;

#define LOG_NUM_THREADS 10
#define NUM_THREADS (1 << LOG_NUM_THREADS)


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
    static const uint32_t TPI=3;                   // threads per instance
    static const uint32_t BITS=512;                 // instance size
    static const uint32_t num_of_bytes=64;                 // instance size

};

typedef cgbn_mem_t<fft_params_t::BITS> Scalar;
typedef cgbn_context_t<fft_params_t::TPI>   context_t;
typedef cgbn_env_t<context_t, fft_params_t::BITS>   env_t;
typedef typename env_t::cgbn_t                bn_t;
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
    const size_t block_length =( 1ul <<  LOG_NUM_THREADS) / fft_params_t::TPI; //TODO lianke when log_m is smaller than log_num_threads,  there is a bug.
    const size_t startidx = idx * block_length;
    if(startidx > length)
        return;
    context_t _context;
    env_t    _env(_context);
    //printf("CGBN with idx=%d, block_length=%d, startidx=%d\n", idx, block_length, startidx);
    env_t::cgbn_t  a, b; 
    /* swapping in place (from Storer's book) */
    for (size_t k = 0; k < block_length; ++k)
    {
        size_t global_k = startidx + k;
        size_t rk = bitreverse(global_k, log_m);
                if (global_k < rk  && rk < length)
        {
            cgbn_load(_env,a, &(input_field[global_k]));
            cgbn_load(_env,b, &(input_field[rk]));
            cgbn_store(_env,( &(input_field[global_k])), b);
            cgbn_store(_env,( &(input_field[rk])), a);
            //cgbn_swap(bn_env, &input_field[global_k], &input_field[rk]);
            // Scalar tmp = input_field[global_k];
            // input_field[global_k] = input_field[rk];
            // input_field[rk] = tmp;
        }
    }
    __syncthreads();
    
}
/*
CGBN API referece. 


void cgbn_modular_power(cgbn_env_t env, cgbn_t &r, const cgbn_t &x, const cgbn_t &e, const cgbn_t &m)
Computes r = x^e modulo the modulus, m. Requires that x < m.


int32_t cgbn_add(cgbn_env_t env, cgbn_t &r, const cgbn_t &a, const cgbn_t &b)
Computes a + b and stores the result in r.   If the sum resulted in a carry out, then 1 is returned to all threads in the group, otherwise return 0.

int32_t cgbn_sub(cgbn_env_t env, cgbn_t &r, const cgbn_t &a, const cgbn_t &b)
Computes a - b and stores the result in r.   If b>a then -1 is returned to all threads, otherwise return 0.


void cgbn_mul(cgbn_env_t env, cgbn_t &r, const cgbn_t &a, const cgbn_t &b)
Computes the low half of the product of a * b, the upper half of the product is discarded.   This is the CGBN equivalent of unsigned multiplication in C.

void cgbn_rem(cgbn_env_t env, cgbn_t &r, const cgbn_t &num, const cgbn_t &denom)
Computes the remainder of num divided by denom and store the result into r, where 0 <= r < denom.
*/

__global__ void cuda_fft_second_step(Scalar *input_field, Scalar omega, const size_t length, const size_t log_m) {
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/fft_params_t::TPI;
    //printf("blockIdx=%d, blockDim.x=%d, threadIdx.x=%d, idx=%d\n",blockIdx.x, blockDim.x, threadIdx.x, idx);
    const size_t block_length =( 1ul <<  LOG_NUM_THREADS) / fft_params_t::TPI; //TODO lianke when log_m is smaller than log_num_threads,  there is a bug.
    const size_t startidx = idx * block_length;
    if(startidx < length){
        size_t m = 1; // invariant: m = 2^{s-1}
        for (size_t s = 1; s <= 1; ++s)
        {
            // w_m is 2^s-th root of unity now
            //TODO lianke need to deal with the data layout difference between CGBN and java.biginteger.tobytearray()
            //TODO lianke replace these arithmetic operators with CGBN.
            const Scalar w_m = omega^(length/(2*m));
            
            for (size_t k = 0; k < block_length; k += 2*m)
            {
                size_t global_k = startidx + k;
                Scalar w = Scalar::one();
                for (size_t j = 0; j < m; ++j)
                {
                    Scalar t = w;
                    t = w * input_field[global_k+j+m];
                    input_field[global_k+j+m] = input_field[global_k+j] - t;
                    input_field[global_k+j] = input_field[global_k+j] + t;
                    w = w * w_m;
                }
            }
            m = m * 2;
        }
    }
    __syncthreads();

}




void best_fft (std::vector<Scalar> &a, const Scalar &omg)
{
	int cnt;
    cudaGetDeviceCount(&cnt);

    printf("CUDA Devices: %d, input_field size: %lu, input_field count: %lu\n", cnt, sizeof(Scalar), a.size());

    size_t threads = NUM_THREADS > 128 ? 128 : NUM_THREADS;
    size_t instance_per_block = (threads / fft_params_t::TPI);//TPI threads per instance, each block has threads.
    size_t blocks = (a.size() + instance_per_block - 1) / instance_per_block;

    printf("NUM_THREADS %u, blocks %lu, threads %lu \n",NUM_THREADS, blocks, threads);
    CUDA_CALL(cudaSetDevice(0));
    Scalar *in; 
    CUDA_CALL( cudaMalloc((void**)&in, sizeof(Scalar) * a.size()); )
    CUDA_CALL( cudaMemcpy(in, (void**)&a[0], sizeof(Scalar) * a.size(), cudaMemcpyHostToDevice); )
    // create a cgbn_error_report for CGBN to report back errors

    const size_t length = a.size();
    const size_t log_m = log2(length); 
    //auto start = std::chrono::steady_clock::now();
    printf("launch block = %d thread = %d\n", blocks, threads);
    cuda_fft_first_step <<<blocks,threads>>>( in, omg, length, log_m);
    CUDA_CALL(cudaDeviceSynchronize());
    cuda_fft_second_step <<<blocks,threads>>>( in, omg, length, log_m);
    CUDA_CALL(cudaDeviceSynchronize());

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

        memcpy(tmp, 
                bytes,
                len);
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

