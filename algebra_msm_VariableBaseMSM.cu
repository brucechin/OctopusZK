#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <tuple>
#include <cmath>
#include "algebra_msm_VariableBaseMSM.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <assert.h>
#include <iostream>
#include <gmp.h>
#include "cgbn/cgbn.h"

#include <bitset>

//TODO lianke: G1 and G2 MSM window table generation can be moved to cpp side too.
using namespace std;

#define CUDA_CALL( call )               \
{                                       \
cudaError_t result = call;              \
if ( cudaSuccess != result )            \
    std::cerr << "CUDA error " << result << " in " << __FILE__ << ":" << __LINE__ << ": " << cudaGetErrorString( result ) << " (" << #call << ")" << std::endl;  \
}


#define CHECK_BIT(var, pos) (((var) >> (pos)) & 1)

class MSM_params_t {
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

//BN254G1_modulus = "21888242871839275222246405745257275088696311157297823662689037894645226208583";
// 811880050|3778125865|3092268470|2172737629|2541841041|1752287885|1008765974|3632069959|
//Fp2Parameters.nonresidue() = "21888242871839275222246405745257275088696311157297823662689037894645226208582";
typedef cgbn_mem_t<MSM_params_t::BITS> Scalar;
typedef cgbn_context_t<MSM_params_t::TPI>   context_t;
typedef cgbn_env_t<context_t, MSM_params_t::BITS>   env_t;
typedef typename env_t::cgbn_t                 bn_t;
typedef typename env_t::cgbn_local_t          bn_local_t;

__device__   uint32_t modulus_raw_G1[16] = {3632069959,1008765974,1752287885,2541841041, 
                                2172737629,3092268470,3778125865,811880050,
                                0, 0,0,0,
                                0,0,0,0};

__device__   uint32_t Fp2_nonresidue_raw[16] = {3632069958,1008765974,1752287885,2541841041, 
                            2172737629,3092268470,3778125865,811880050,
                            0, 0,0,0,
                            0,0,0,0};   

__device__   uint32_t zero_raw[16] = {0,0,0,0, 
                                      0,0,0,0,
                                      0, 0,0,0,
                                      0,0,0,0};
// Declare the instance type
typedef struct {
  cgbn_mem_t<MSM_params_t::BITS> X;
  cgbn_mem_t<MSM_params_t::BITS> Y;
  cgbn_mem_t<MSM_params_t::BITS> Z;
} BN254G1;//this is raw memory struct.

typedef struct {
  cgbn_mem_t<MSM_params_t::BITS> Xa;
  cgbn_mem_t<MSM_params_t::BITS> Xb;
  cgbn_mem_t<MSM_params_t::BITS> Ya;
  cgbn_mem_t<MSM_params_t::BITS> Yb;
  cgbn_mem_t<MSM_params_t::BITS> Za;
  cgbn_mem_t<MSM_params_t::BITS> Zb;
} BN254G2;

typedef struct {
  bn_t X;
  bn_t Y;
  bn_t Z;
} BN254G1Compute;

typedef struct {
  bn_t a;
  bn_t b;
} Fp2;

typedef struct {
  Fp2 X;
  Fp2 Y;
  Fp2 Z;
} BN254G2Compute;

__device__ __forceinline__
bool testBit(Scalar input, int n)
{
    int byte_index =n / 32;
    int byte_offset = n % 32;
    return CHECK_BIT(input._limbs[byte_index], byte_offset);
}



__device__ __forceinline__
bool isZero(BN254G1Compute input)
{
    context_t _context;
    env_t    _env(_context);
    bn_t zero;
    Scalar zero_binary;
    memcpy(zero_binary._limbs, zero_raw, MSM_params_t::num_of_bytes);
    cgbn_load(_env, zero, &zero_binary);
    
    //printf("input last uint=%d, isZero=%d\n", input.Z._limbs[0], cgbn_equals(_env, zero, z));
    return cgbn_equals(_env, zero, input.Z);
}

__device__ __forceinline__
bool isZero(BN254G2Compute input)
{
    context_t _context;
    env_t    _env(_context);
    bn_t zero;
    Scalar zero_binary;
    memcpy(zero_binary._limbs, zero_raw, MSM_params_t::num_of_bytes);
    cgbn_load(_env, zero, &zero_binary);
    
    //printf("input last uint=%d, isZero=%d\n", input.Z._limbs[0], cgbn_equals(_env, zero, z));
    return cgbn_equals(_env, zero, input.Z.a) && cgbn_equals(_env, zero, input.Z.b);
}

__device__ __forceinline__
BN254G1Compute twice(BN254G1Compute a)
{
  context_t _context;
  env_t    _env(_context);
  if(isZero(a)){
    return a;
  }

  Scalar modulus_binary;
  memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
  bn_t modulus;
  cgbn_load(_env, modulus, &modulus_binary);

  BN254G1Compute result;
  memset(result.X._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Y._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Z._limbs, 0, MSM_params_t::num_of_bytes);
  bn_t a_x, a_y, a_z;
  a_x = a.X;
  a_y = a.Y;
  a_z = a.Z;

  bn_t A,B,C,D,E,F,X3,Y3,Y1Z1,Z3, eightC;


  cgbn_mul(_env, A, a_x, a_x);
  cgbn_rem(_env, A, A, modulus);

  cgbn_mul(_env, B, a_y, a_y);
  cgbn_rem(_env, B, B, modulus);


  cgbn_mul(_env, C, B, B);
  cgbn_rem(_env, C, C, modulus);


 // D = 2 * ((X1 + B)^2 - A - C)
  cgbn_add(_env, D, a_x, B);
  cgbn_rem(_env, D, D, modulus);
  cgbn_mul(_env, D, D, D);
  cgbn_rem(_env, D, D, modulus);
  cgbn_add(_env, D, D, modulus);
  cgbn_add(_env, D, D, modulus);
  cgbn_sub(_env, D, D, A);
  cgbn_sub(_env, D, D, C);
  cgbn_add(_env, D, D, D);
  cgbn_rem(_env, D, D, modulus);

  // E = 3 * A
  cgbn_add(_env, E, A, A);
  cgbn_add(_env, E, E, A);
  cgbn_rem(_env, E, E, modulus);

  // F = E^2
  cgbn_mul(_env, F, E, E);
  cgbn_rem(_env, F, F, modulus);



   // X3 = F - 2 D
  cgbn_add(_env, X3, F, modulus);
  cgbn_add(_env, X3, F, modulus);
  cgbn_sub(_env, X3, X3, D);
  cgbn_sub(_env, X3, X3, D);
  cgbn_rem(_env, X3, X3, modulus);




  cgbn_add(_env, eightC, C, C);
  cgbn_add(_env, eightC, eightC, eightC);
  cgbn_add(_env, eightC, eightC, eightC);
  cgbn_rem(_env, eightC, eightC, modulus);


  // Y3 = E * (D - X3) - 8 * C
  cgbn_add(_env, Y3, D, modulus);
  cgbn_sub(_env, Y3, Y3, X3);
  cgbn_rem(_env, Y3, Y3, modulus);
  cgbn_mul(_env, Y3, Y3, E);
  cgbn_rem(_env, Y3, Y3, modulus);
  cgbn_add(_env, Y3, Y3, modulus);
  cgbn_sub(_env, Y3, Y3, eightC);
  cgbn_rem(_env, Y3, Y3, modulus);


  // Z3 = 2 * Y1 * Z1
  cgbn_mul(_env, Y1Z1, a_y, a_z);
  cgbn_rem(_env, Y1Z1, Y1Z1, modulus);
  cgbn_add(_env, Z3, Y1Z1, Y1Z1);
  cgbn_rem(_env, Z3, Z3, modulus);

  // cgbn_store(_env, &result.X, X3);
  // cgbn_store(_env, &result.Y, Y3);
  // cgbn_store(_env, &result.Z, Z3);
  result.X = X3;
  result.Y = Y3;
  result.Z = Z3;

  return result;

}



__device__ __forceinline__
BN254G1Compute add(BN254G1Compute a, BN254G1Compute b) {
    // // Handle special cases having to do with O

  // printf("11111");
  if (isZero(a)) {
      return b;
  }

  if (isZero(b)) {
      return a;
  }
  // printf("22222");

  context_t _context;
  env_t    _env(_context);


  Scalar modulus_binary;
  memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
  bn_t modulus;
  cgbn_load(_env, modulus, &modulus_binary);

  BN254G1Compute result;
  memset(result.X._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Y._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Z._limbs, 0, MSM_params_t::num_of_bytes);

  bn_t a_x, a_y, a_z, b_x, b_y, b_z;
  a_x = a.X; 
  a_y = a.Y; 
  a_z = a.Z;
  b_x = b.X;
  b_y = b.Y;
  b_z = b.Z;


  bn_t Z1Z1, Z2Z2, U1, U2, Z1_cubed, Z2_cubed, S1, S2;
  
  

  cgbn_mul(_env, Z1Z1, a_z, a_z);
  cgbn_rem(_env, Z1Z1, Z1Z1, modulus);

  cgbn_mul(_env, Z2Z2, b_z, b_z);
  cgbn_rem(_env, Z2Z2, Z2Z2, modulus);

  cgbn_mul(_env, U1, a_x, Z2Z2); 
  cgbn_rem(_env, U1, U1, modulus);

  cgbn_mul(_env, U2, b_x, Z1Z1);
  cgbn_rem(_env, U2, U2, modulus);

  cgbn_mul(_env, Z1_cubed, a_z, Z1Z1);
  cgbn_rem(_env, Z1_cubed, Z1_cubed, modulus);

  cgbn_mul(_env, Z2_cubed, b_z, Z2Z2);
  cgbn_rem(_env, Z2_cubed, Z2_cubed, modulus);

  cgbn_mul(_env, S1, a_y, Z2_cubed);
  cgbn_rem(_env, S1, S1, modulus);

  cgbn_mul(_env, S2, b_y, Z1_cubed);
  cgbn_rem(_env, S2, S2, modulus);
  // printf("333333");


  if (cgbn_equals(_env, U1, U2) && cgbn_equals(_env, S1, S2)) {
      // Double case; nothing above can be reused.
      return twice(a);
  }
  // printf("444444");

  bn_t H, S2_minus_S1, I, J, r, V, X3, S1_J, Y3, Z3;
  
  // H = U2-U1
  cgbn_add(_env, H, U2, modulus);
  cgbn_sub(_env, H, H, U1);
  cgbn_rem(_env, H, H, modulus);

  cgbn_add(_env, S2_minus_S1, S2, modulus);
  cgbn_sub(_env, S2_minus_S1, S2_minus_S1, S1);
  cgbn_rem(_env, S2_minus_S1, S2_minus_S1, modulus);


  // I = (2 * H)^2
  cgbn_add(_env, I, H, H);
  cgbn_rem(_env, I, I, modulus);
  cgbn_mul(_env, I, I, I);
  cgbn_rem(_env, I, I, modulus);

  // J = H * I
  cgbn_mul(_env, J, H, I);
  cgbn_rem(_env, J, J, modulus);

  // r = 2 * (S2-S1)
  cgbn_add(_env, r, S2_minus_S1, S2_minus_S1);
  cgbn_rem(_env, r, r, modulus);

  // V = U1 * I
  cgbn_mul(_env, V, U1, I);
  cgbn_rem(_env, V, V, modulus);

  // X3 = r^2 - J - 2 * V
  cgbn_mul(_env, X3, r, r);
  cgbn_rem(_env, X3, X3, modulus);
  cgbn_add(_env, X3, X3, modulus);
  cgbn_add(_env, X3, X3, modulus);
  cgbn_add(_env, X3, X3, modulus);
  cgbn_sub(_env, X3, X3, J);
  cgbn_sub(_env, X3, X3, V);
  cgbn_sub(_env, X3, X3, V);
  cgbn_rem(_env, X3, X3, modulus);


  // Y3 = r * (V-X3)-2 * S1_J
  cgbn_mul(_env, S1_J, S1, J);
  cgbn_rem(_env, S1_J, S1_J, modulus);
  cgbn_add(_env, Y3, V, modulus);
  cgbn_sub(_env, Y3, Y3, X3);
  cgbn_rem(_env, Y3, Y3, modulus);
  cgbn_mul(_env, Y3, Y3, r);
  cgbn_rem(_env, Y3, Y3, modulus);
  cgbn_add(_env, Y3, Y3, modulus);
  cgbn_add(_env, Y3, Y3, modulus);
  cgbn_sub(_env, Y3, Y3, S1_J);
  cgbn_sub(_env, Y3, Y3, S1_J);
  cgbn_rem(_env, Y3, Y3, modulus);



  cgbn_add(_env, Z3, a_z, b_z);
  cgbn_rem(_env, Z3, Z3, modulus);
  cgbn_mul(_env, Z3, Z3, Z3);
  cgbn_rem(_env, Z3, Z3, modulus);
  cgbn_add(_env, Z3, Z3, modulus);
  cgbn_add(_env, Z3, Z3, modulus);
  cgbn_sub(_env, Z3, Z3, Z1Z1);
  cgbn_sub(_env, Z3, Z3, Z2Z2);
  cgbn_rem(_env, Z3, Z3, modulus);
  cgbn_mul(_env, Z3, Z3, H);
  cgbn_rem(_env, Z3, Z3, modulus);
  // printf("555555\n");

  result.X = X3;
  result.Y = Y3;
  result.Z = Z3;

  return result;

}


__global__ void pippengerMSMG1_unit(Scalar* inputScalarArray, BN254G1* inputBaseArray, BN254G1* outputBN254Array, int batch_size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    context_t _context;
    env_t    _env(_context);
    BN254G1Compute res;
    cgbn_load(_env, res.X, &inputBaseArray[0].X);
    cgbn_load(_env, res.Y, &inputBaseArray[0].Y);
    cgbn_load(_env, res.Z, &inputBaseArray[0].Z);

    if(idx >= batch_size){
      return;
    }


    //TODO it seems that pippengerMSM is not easy to parallelize.

    // for (int outer = 0; outer < outerc; ++outer) {
    //     int inner = 0;
    //     for (int i = 0; i < windowSize; ++i) {
    //         //testBit is correct
    //         if (testBit(inputScalarArray[idx], outer * windowSize + i)) {
    //             inner |= 1 << i;
    //         }
    //     }
    //     //printf("idx=%d, outer=%d, inner=%d\n", idx, outer, inner);
    //     BN254G1Compute to_added;
    //     cgbn_load(_env, to_added.X, &inputBaseArray[tableInnerSize * outer + inner].X);
    //     cgbn_load(_env, to_added.Y, &inputBaseArray[tableInnerSize * outer + inner].Y);
    //     cgbn_load(_env, to_added.Z, &inputBaseArray[tableInnerSize * outer + inner].Z);
    //     res = add(res, to_added);
    // }



    cgbn_store(_env, &outputBN254Array[idx].X, res.X);
    cgbn_store(_env, &outputBN254Array[idx].Y, res.Y);
    cgbn_store(_env, &outputBN254Array[idx].Z, res.Z);

    return;
}

void  pippengerMSMG1(std::vector<Scalar> & bigScalarArray, std::vector<BN254G1> &multiplesOfBasePtrArray, BN254G1* outputArray)
{
	int cnt;
    cudaGetDeviceCount(&cnt);
    size_t batch_size = bigScalarArray.size();

    printf("CUDA Devices: %d, input_field size: %lu, input_field count: %lu\n", cnt, sizeof(Scalar), batch_size);
    size_t threads_per_block = 128;
    size_t instance_per_block = (threads_per_block / MSM_params_t::TPI);//TPI threads per instance, each block has threads.
    size_t blocks = (batch_size + instance_per_block - 1) / instance_per_block;
    printf("num of blocks %lu, threads per block %lu \n", blocks, threads_per_block);
    CUDA_CALL(cudaSetDevice(0));
    Scalar *inputScalarArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputScalarArrayGPU, sizeof(Scalar) * batch_size); )
    CUDA_CALL( cudaMemcpy(inputScalarArrayGPU, (void**)&bigScalarArray[0], sizeof(Scalar) * batch_size, cudaMemcpyHostToDevice); )
    
    BN254G1 *inputBaseArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputBaseArrayGPU, sizeof(BN254G1) * batch_size); )
    CUDA_CALL( cudaMemcpy(inputBaseArrayGPU, (void**)&multiplesOfBasePtrArray[0], sizeof(BN254G1) * batch_size, cudaMemcpyHostToDevice); )
    
    BN254G1* outputBN254ArrayGPU;//TODO this is partial results, should be aggregated again to obtain one G1 value.
    CUDA_CALL( cudaMalloc((void**)&outputBN254ArrayGPU, sizeof(BN254G1) * batch_size); )
    CUDA_CALL( cudaMemset(outputBN254ArrayGPU, 0, sizeof(BN254G1) * batch_size); )

    printf("launch block = %d thread = %d\n", blocks, threads_per_block);


    int numBits = 254;//BN254 specific value;
    int length = batch_size;
    int log2Length = max(1, log2(length));
    int c= log2Length - (log2Length/3);
    int numBuckets = 1 << c;
    int numGroups = (numBits + c - 1 ) / c;
    BN254G1  zero;
    zero.Y._limbs[0] = 1;
    vector<BN254G1> buketsModel(numBuckets, zero);
    //TODO lianke
    for(int k = numGroups - 1; k >= 0; k--){

        pippengerMSMG1_unit <<<blocks,threads_per_block, 32 * 1024>>>( inputScalarArrayGPU, inputBaseArrayGPU, outputBN254ArrayGPU, batch_size);
        CUDA_CALL(cudaDeviceSynchronize());
    }

    //TODO lianke write another cuda kernel to aggregate outputBN254ArrayGPU.


    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }

    CUDA_CALL(cudaMemcpy((void**)outputArray, outputBN254ArrayGPU, sizeof(BN254G1) * batch_size, cudaMemcpyDeviceToHost); )
    CUDA_CALL(cudaDeviceSynchronize());

    CUDA_CALL(cudaFree(inputScalarArrayGPU));
    CUDA_CALL(cudaFree(inputBaseArrayGPU));
    CUDA_CALL(cudaFree(outputBN254ArrayGPU));

}




/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseSerialMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseSerialMSMNativeHelper
  (JNIEnv * env, jclass obj, jobject multiplesOfBaseX, jobject multiplesOfBaseY, jobject multiplesOfBaseZ, jobject scalars){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");


    jint base_size = env->CallIntMethod(bases, java_util_ArrayList_size);
    jint scalars_size = env->CallIntMethod(scalars, java_util_ArrayList_size);


    vector<Scalar> bigScalarArray = vector<Scalar>(batch_size, Scalar());
    vector<BN254G1> multiplesOfBasePtrArray = vector<BN254G1>(base_size, BN254G1());


    for(int i =0; i < batch_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)&bigScalarArray[i]._limbs;
        memcpy(tmp, bytes, len);
    }

    for(int i = 0; i < base_size;i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(multiplesOfBaseX, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)multiplesOfBasePtrArray[i].X._limbs;
        memcpy(tmp, bytes,len);

        
    }

    for(int i = 0; i < base_size;i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(multiplesOfBaseY, java_util_ArrayList_get, i) ;
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)multiplesOfBasePtrArray[i].Y._limbs;
        memcpy(tmp, bytes,len);
        
    }

    for(int i = 0; i < base_size;i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(multiplesOfBaseZ, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        int len = env->GetArrayLength(element);
        char* tmp = (char*)multiplesOfBasePtrArray[i].Z._limbs;
        memcpy(tmp, bytes,len);
        
    }


    //lianke: we know in advance that the numBits will be at most 254. we hard encode it.



    jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes);


    int length = filteredInput.size();
    int log2Length =  max(1, (int)log2(length));
    int c = log2Length - (log2Length / 3);
    int numBuckets = 1 << c;
    int numGroups = (numBits + c - 1)/c;
    BigInt zero("0"); 
    vector<BigInt> bucketsModel = vector<BigInt>(numBuckets, zero);
    BigInt result("0"); 
    for(int k = numGroups - 1; k >=0; k--){
        if (k < numGroups - 1) {
            for (int i = 0; i < c; i++) {
                result = (result + result) ;
                result %= FqModulusParameter; 
            }
        }
        vector<BigInt> buckets = vector<BigInt>(bucketsModel);

    for (int i = 0; i < length; i++) {
            int id = 0;
            for (int j = 0; j < c; j++) {
                if (std::get<0>(filteredInput[i]).testBit(k * c + j)) {
                    id |= 1 << j;
                }
            }

            if (id == 0) {
                continue;
            }

            // Potentially use mixed addition here.
            buckets[id] = (buckets[id] + std::get<1>(filteredInput[i]));
            buckets[id]  %= FqModulusParameter;
    }


    BigInt runningSum("0");
    for(int i = numBuckets - 1; i > 0; i--){
        runningSum = (runningSum + buckets[i]);
        runningSum %= FqModulusParameter;
        result = (result + runningSum);
        result %= FqModulusParameter;
    }
    }


    acc = (acc + result);
    acc %= FqModulusParameter;

    env->SetByteArrayRegion(resultByteArray, 0 , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(acc.bytes));

    return resultByteArray;

  }