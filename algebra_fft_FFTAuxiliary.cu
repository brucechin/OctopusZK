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
    static const uint32_t TPI=32;                   // threads per instance
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
    uint32_t Fr_modulus_raw[16] = {4026531841,1138881939,2042196113,674490440, 
                                2172737629,3092268470,3778125865,811880050,
                                0, 0,0,0,
                                0,0,0,0};
    
    
    memcpy(modulus_binary._limbs, Fr_modulus_raw, fft_params_t::num_of_bytes);
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

}




void best_fft (std::vector<Scalar> &a, const Scalar &omg, int taskID)
{
    int num_gpus = 1;
    CUDA_CALL(cudaGetDeviceCount(&num_gpus));

    size_t threads_per_block = 128;
    size_t instance_per_block = (threads_per_block / fft_params_t::TPI);//TPI threads per instance, each block has threads.
    size_t blocks = (a.size() + instance_per_block - 1) / instance_per_block;

    CUDA_CALL(cudaSetDevice(taskID % num_gpus));
    Scalar *in; 
    CUDA_CALL( cudaMalloc((void**)&in, sizeof(Scalar) * a.size()); )
    CUDA_CALL( cudaMemcpy(in, (void**)&a[0], sizeof(Scalar) * a.size(), cudaMemcpyHostToDevice); )

    const size_t length = a.size();
    const size_t log_m = log2(length); 

    cuda_fft_first_step <<<blocks,threads_per_block>>>( in, omg, length, log_m);
    CUDA_CALL(cudaDeviceSynchronize());
    size_t s = 1;
    for(; s <= log_m; s++){
        cuda_fft_second_step <<<blocks,threads_per_block>>>( in, omg, length, log_m, s);
        CUDA_CALL(cudaDeviceSynchronize());
    }

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


__global__ void cuda_fft_first_step_batch( Scalar *input_field, Scalar omega, const size_t batch_size, const size_t instance_size, const size_t log_m) {

    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/fft_params_t::TPI;
    //printf("blockIdx=%d, blockDim.x=%d, threadIdx.x=%d, idx=%d\n",blockIdx.x, blockDim.x, threadIdx.x, idx);
    //const size_t block_length =( 1ul <<  LOG_NUM_THREADS) / fft_params_t::TPI;
    if(idx > batch_size * instance_size)
        return;
    context_t _context;
    env_t    _env(_context);

    bn_t  a, b; 
    /* swapping in place (from Storer's book) */
    size_t group_id = idx / instance_size;
    size_t k_within_instance = idx % instance_size;
    size_t rk_within_instance = bitreverse(k_within_instance, log_m);
    size_t global_rk = group_id * instance_size + rk_within_instance;
    size_t global_k = idx;
    if (global_k < global_rk  && rk_within_instance < instance_size)
    {
        cgbn_load(_env,a, &(input_field[global_k]));
        cgbn_load(_env,b, &(input_field[global_rk]));
        cgbn_store(_env, &(input_field[global_k]), b);
        cgbn_store(_env, &(input_field[global_rk]), a);
    }
    
}


__global__ void cuda_fft_second_step_batch(Scalar *input_field, Scalar omega_binary, const size_t batch_size, const size_t instance_size, const size_t log_m, size_t s_index) {
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/fft_params_t::TPI;

    context_t _context;
    env_t    _env(_context);

    Scalar modulus_binary;
    //Fr modulus is : |811880050|3778125865|3092268470|2172737629|674490440|2042196113|1138881939|4026531841|
    uint32_t Fr_modulus_raw[16] = {4026531841,1138881939,2042196113,674490440, 
                                2172737629,3092268470,3778125865,811880050,
                                0, 0,0,0,
                                0,0,0,0};
    
    
    memcpy(modulus_binary._limbs, Fr_modulus_raw, fft_params_t::num_of_bytes);
    bn_t modulus;
    cgbn_load(_env, modulus, &modulus_binary);
    bn_t omega;
    cgbn_load(_env, omega, &omega_binary);

    size_t m = 1 << (s_index - 1); // invariant: m = 2^{s-1}
    Scalar exponential_binary;
    exponential_binary._limbs[0] = (uint32_t)instance_size/(2*m);
    bn_t exponential, w_m;
    cgbn_load(_env, exponential, &exponential_binary);
    cgbn_modular_power(_env, w_m, omega, exponential, modulus);
    // w_m is 2^s-th root of unity now


    size_t group_id = idx / instance_size;
    size_t k_within_instance = idx % instance_size;
    size_t k_plus_j_within_instance = (k_within_instance / m) * m * 2 + k_within_instance %  m;

    size_t global_kj = group_id * instance_size + k_plus_j_within_instance;
    if(k_plus_j_within_instance < instance_size){
        bn_t w, w_exp;
        size_t w_exp_int = idx % m;
        Scalar w_exp_binary;
        w_exp_binary._limbs[0] = w_exp_int;
        cgbn_load(_env, w_exp, &w_exp_binary);
        cgbn_modular_power(_env, w, w_m, w_exp, modulus);

        bn_t t;
        bn_t input_kjm, input_kj;
        cgbn_load(_env, input_kjm, &input_field[global_kj +m]);
        cgbn_load(_env, input_kj, &input_field[global_kj]);
        cgbn_mul(_env, t, w, input_kjm);
        cgbn_rem(_env, t, t, modulus);

        cgbn_add(_env, input_kj, input_kj, modulus);
        cgbn_sub(_env, input_kj, input_kj, t);
        //after subtraction, the result could be negative. so i need to add modulus to make it a positive number 
        cgbn_rem(_env, input_kj, input_kj, modulus);

        cgbn_store(_env, &input_field[global_kj + m], input_kj);
        cgbn_load(_env, input_kj, &input_field[global_kj]);
        cgbn_add(_env, input_kj, input_kj, t);
        cgbn_rem(_env, input_kj, input_kj, modulus);

        cgbn_store(_env, &input_field[global_kj], input_kj);
        cgbn_mul(_env, w, w, w_m);
        cgbn_rem(_env, w, w, modulus);
    }

}


//we are process batch_size small FFT tasks, each task contains instance_size elements.
void best_fft_batch (std::vector<Scalar> &a, const Scalar &omg, int batch_size, int instance_size, int taskID)
{
    int num_gpus = 1;
    CUDA_CALL(cudaGetDeviceCount(&num_gpus));
    size_t threads_per_block = 128;
    size_t instance_per_block = (threads_per_block / fft_params_t::TPI);//TPI threads per instance, each block has threads.
    size_t blocks = (a.size() + instance_per_block - 1) / instance_per_block;

    CUDA_CALL(cudaSetDevice(taskID % num_gpus));
    Scalar *in; 
    CUDA_CALL( cudaMalloc((void**)&in, sizeof(Scalar) * a.size()); )
    CUDA_CALL( cudaMemcpy(in, (void**)&a[0], sizeof(Scalar) * a.size(), cudaMemcpyHostToDevice); )

    const size_t log_m = log2(instance_size); 

    cuda_fft_first_step_batch <<<blocks,threads_per_block>>>( in, omg, batch_size, instance_size, log_m);
    CUDA_CALL(cudaDeviceSynchronize());
    size_t s = 1;
    for(; s <= log_m; s++){
        cuda_fft_second_step_batch <<<blocks,threads_per_block>>>( in, omg, batch_size, instance_size, log_m, s);
        CUDA_CALL(cudaDeviceSynchronize());
    }

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
  (JNIEnv * env, jclass obj, jobject inputs, jbyteArray omegaArray, jint taskID){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");
    jint input_len = env->CallIntMethod(inputs, java_util_ArrayList_size);

    vector<Scalar> inputArray = vector<Scalar>(input_len, Scalar());
    for(int i =0; i < input_len; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(inputs, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)&inputArray[i]._limbs;

        memcpy(tmp, bytes, len);
    }

    //cout << "obtain all input arrays" << endl;
    Scalar omega;
    char* bytes = (char*)env->GetByteArrayElements(omegaArray, NULL);
    int len = env->GetArrayLength(omegaArray);
    char* tmp = (char*)&omega._limbs;
    memcpy(tmp , bytes, len);
    
    best_fft(inputArray, omega, taskID);

    //cout << "finish cpp fft on input arrays" << endl;
    //TODO lianke when the size is too large, this will over flow. we should do ArrayList<byte[]> although it may be slower.
    jbyteArray resultByteArray = env->NewByteArray((jsize)fft_params_t::num_of_bytes * input_len);
    env->SetByteArrayRegion(resultByteArray, 0 , fft_params_t::num_of_bytes * input_len,   reinterpret_cast< jbyte*>(&inputArray[0]));
    
    return resultByteArray;

}

/*
 * Class:     algebra_fft_FFTAuxiliary
 * Method:    serialRadix2FFTNativeBatchPartition
 * Signature: (Ljava/util/List;[BIII)Ljava/util/ArrayList;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_fft_FFTAuxiliary_serialRadix2FFTNativeBatchPartition
  (JNIEnv * env, jclass obj, jobject inputs, jbyteArray omegaArray, jint batch_size, jint instance_size, jint taskID){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");

    int len = 32;
    vector<Scalar> inputArray = vector<Scalar>(batch_size * instance_size, Scalar());
    for(int i =0; i < batch_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(inputs, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        for(int j = 0; j < instance_size; j++){
            char* tmp = (char*)&inputArray[i * instance_size + j];
            memcpy(tmp, &bytes[j * len], len);
        }
    }

    Scalar omega;
    char* bytes = (char*)env->GetByteArrayElements(omegaArray, NULL);
    char* tmp = (char*)&omega._limbs;
    memcpy(tmp , bytes, len);

    best_fft_batch(inputArray, omega, batch_size, instance_size, taskID);
   //cout << "CUDA FFT finished" << endl;

    jbyteArray resultByteArray = env->NewByteArray((jsize) fft_params_t::num_of_bytes * instance_size * batch_size);
    env->SetByteArrayRegion(resultByteArray, 0 , fft_params_t::num_of_bytes * instance_size * batch_size,   reinterpret_cast< jbyte*>(&inputArray[0]));



    return resultByteArray;

  }
