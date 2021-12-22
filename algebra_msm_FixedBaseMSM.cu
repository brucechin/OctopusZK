#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <chrono>
#include "algebra_msm_FixedBaseMSM.h"
//#include "BigInteger.h"
//#include "BigInt.h"
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <assert.h>
#include <vector>
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
typedef cgbn_mem_t<MSM_params_t::BITS> Scalar;
typedef cgbn_context_t<MSM_params_t::TPI>   context_t;
typedef cgbn_env_t<context_t, MSM_params_t::BITS>   env_t;
typedef typename env_t::cgbn_t                 bn_t;
typedef typename env_t::cgbn_local_t          bn_local_t;

__device__   uint32_t modulus_raw_G1[16] = {3632069959,1008765974,1752287885,2541841041, 
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
} BN254G1;

// typedef struct {
//   cgbn_mem_t<BITS> X;
//   cgbn_mem_t<BITS> Y;
//   cgbn_mem_t<BITS> Z;
// } BN254G2;



__device__ __forceinline__
bool testBit(Scalar input, int n)
{
    int byte_index =n / 32;
    int byte_offset = n % 32;
    return CHECK_BIT(input._limbs[byte_index], byte_offset);
}


__device__ __forceinline__
bool isZero(BN254G1 input)
{
    context_t _context;
    env_t    _env(_context);
    bn_t zero, z;
    Scalar zero_binary;
    memcpy(zero_binary._limbs, zero_raw, MSM_params_t::num_of_bytes);
    cgbn_load(_env, zero, &zero_binary);
    cgbn_load(_env, z, &input.Z);
    
    //printf("input last uint=%d, isZero=%d\n", input.Z._limbs[0], cgbn_equals(_env, zero, z));
    return cgbn_equals(_env, zero, z);
}

__device__ __forceinline__
BN254G1 twice(BN254G1 a)
{
  context_t _context;
  env_t    _env(_context);


  Scalar modulus_binary;
  memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
  bn_t modulus;
  cgbn_load(_env, modulus, &modulus_binary);

  BN254G1 result;
  bn_t a_x, a_y, a_z;
  cgbn_load(_env, a_x, &a.X);
  cgbn_load(_env, a_y, &a.Y);
  cgbn_load(_env, a_z, &a.Z);

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
  cgbn_add(_env, E, A, A);
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

  cgbn_store(_env, &result.X, X3);
  cgbn_store(_env, &result.Y, Y3);
  cgbn_store(_env, &result.Z, Z3);

  return result;

}


__device__ __forceinline__
bool equals(BN254G1 a, BN254G1 b)
{
  if(isZero(a)){
    return isZero(b);
  }

  if(isZero(b)){
    return false;
  }

  context_t _context;
  env_t    _env(_context);

  Scalar modulus_binary;
  memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
  bn_t modulus;
  cgbn_load(_env, modulus, &modulus_binary);

  bn_t a_x, a_y, a_z, b_x, b_y, b_z;
  cgbn_load(_env, a_x, &a.X);
  cgbn_load(_env, a_y, &a.Y);
  cgbn_load(_env, a_z, &a.Z);
  cgbn_load(_env, b_x, &b.X);
  cgbn_load(_env, b_y, &b.Y);
  cgbn_load(_env, b_z, &b.Z);   


  bn_t Z1_squared, Z2_squared, XZ1_squared, XZ2_squared, Z1_cubed, Z2_cubed;
  cgbn_mul(_env, Z1_squared, a_z, a_z);
  cgbn_rem(_env, Z1_squared, Z1_squared, modulus);

  cgbn_mul(_env, Z2_squared, b_z, b_z);
  cgbn_rem(_env, Z2_squared, Z2_squared, modulus);

  cgbn_mul(_env, XZ1_squared, Z1_squared, b_x);
  cgbn_rem(_env, XZ1_squared, XZ1_squared, modulus);

  cgbn_mul(_env, XZ2_squared, Z2_squared, a_x);
  cgbn_rem(_env, XZ2_squared, XZ2_squared, modulus);


  if(cgbn_equals(_env, XZ1_squared, XZ2_squared)){
    return false;
  }

  cgbn_mul(_env, Z1_cubed, Z1_squared, a_z);
  cgbn_rem(_env, Z1_cubed, Z1_cubed, modulus);

  cgbn_mul(_env, Z2_cubed, Z2_squared, b_z);
  cgbn_rem(_env, Z2_cubed, Z2_cubed, modulus);


  bn_t YZ2_cubed, YZ1_cubed;
  cgbn_mul(_env, YZ1_cubed, Z1_cubed, b_y);
  cgbn_rem(_env, YZ1_cubed, YZ1_cubed, modulus);

  cgbn_mul(_env, YZ2_cubed, Z2_cubed, a_y);
  cgbn_rem(_env, YZ2_cubed, YZ2_cubed, modulus);

  if(cgbn_equals(_env, YZ1_cubed, YZ2_cubed)){
    return false;
  }

  return true;
}

__device__ __forceinline__
void printMem(Scalar input)
{
    for(int i = 0; i < MSM_params_t::BITS/32; i++){
      printf("%lu|", input._limbs[i]);
    }
    printf("finished\n");
}


__device__ __forceinline__
void print_bn_t(bn_t &number) {
  using __env_t = bn_t::parent_env_t;
  const int IPB = blockDim.x/__env_t::TPI
  __shared__ uint32_t n[IPB][(__env_t::BITS/32)] ;
  __shared__ uint32_t vote[IPB];
  bool is_represent = (threadIdx.x % TPI) == 0;
  bool instance_id  = threadIdx.x / TPI;
  bool tid_in_instance = threadIdx.x % TPI;
  if (is_represent) vote[instance_id] = 0;
  for (int i = 0; i < __env_t::LIMBS; i++)
    n[instance_id][tid_in_instance * __env_t::LIMBS + i] = number._limbs[i];
  atomicAdd(&vote[instance_id], 1);
  while (vote[instance_id] < TPI) ;
  if (is_represent) {
    printf("instance %d is ", (threadIdx.x + blockIdx.x * blockDim.x)/TPI );
    for (int i = 0; i < __env_t::BITS/32; i++) {
      printf(" %X |", n[instance_id][i]);
    }
    printf("\n");
  }
}

__device__ __forceinline__
BN254G1 add(BN254G1 a, BN254G1 b) {
    // // Handle special cases having to do with O

    
  if (isZero(a)) {
      return b;
  }

  if (isZero(b)) {
      return a;
  }
  context_t _context;
  env_t    _env(_context);


  Scalar modulus_binary;
  memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
  bn_t modulus;
  cgbn_load(_env, modulus, &modulus_binary);

  BN254G1 result;
  bn_t a_x, a_y, a_z, b_x, b_y, b_z;
  cgbn_load(_env, a_x, &a.X);
  cgbn_load(_env, a_y, &a.Y);
  cgbn_load(_env, a_z, &a.Z);
  cgbn_load(_env, b_x, &b.X);
  cgbn_load(_env, b_y, &b.Y);
  cgbn_load(_env, b_z, &b.Z);

  bn_t Z1Z1, Z2Z2, U1, U2, Z1_cubed, Z2_cubed, S1, S2;
  
  
  //Scalar tmp;
  //cgbn_store(_env, &tmp, a_z);
  //printf("z=");
  //printMem(tmp);
  cgbn_mul(_env, Z1Z1, a_z, a_z);
  //cgbn_store(_env, &tmp, Z1Z1);
  //printf("Z1Z1=");
  //printMem(tmp);
  cgbn_rem(_env, Z1Z1, Z1Z1, modulus);
  //cgbn_store(_env, &tmp, Z1Z1);
  //printf("Z1Z1 mod=");
  //printMem(tmp);
  // printf("\n\n");
  // printf("\n\n");
  // printf("\n\n");


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


  if (cgbn_equals(_env, U1, U2) && cgbn_equals(_env, S1, S2)) {
      // Double case; nothing above can be reused.
      return twice(a);
  }

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

  cgbn_store(_env, &result.X, X3);
  cgbn_store(_env, &result.Y, Y3);
  cgbn_store(_env, &result.Z, Z3);

  return result;

}

__global__ void MSM_unit_processing(Scalar* inputScalarArray, BN254G1* inputBaseArray, BN254G1* outputBN254Array, int outerc, int windowSize, int tableInnerSize, int batch_size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    BN254G1 res = inputBaseArray[0];
    if(idx >= batch_size){
      //printf("found large idx=%d\n", idx);
      return;
    }
    for (int outer = 0; outer < outerc; ++outer) {
        int inner = 0;
        for (int i = 0; i < windowSize; ++i) {
            //testBit is correct
            if (testBit(inputScalarArray[idx], outer * windowSize + i)) {
                inner |= 1 << i;
            }
        }
        //printf("CUDA outer=%d, inner=%d\n", outer, inner);
        //TODO add is wrong
        res = add(res, inputBaseArray[tableInnerSize * outer + inner]);
    }
    outputBN254Array[idx] = res;

}


void  fixed_batch_MSM(std::vector<Scalar> & bigScalarArray, std::vector<BN254G1> &multiplesOfBasePtrArray, BN254G1* outputArray,int outerc, int windowSize, int out_len, int inner_len)
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
    CUDA_CALL( cudaMalloc((void**)&inputBaseArrayGPU, sizeof(BN254G1) * out_len * inner_len); )
    CUDA_CALL( cudaMemcpy(inputBaseArrayGPU, (void**)&multiplesOfBasePtrArray[0], sizeof(BN254G1) * out_len * inner_len, cudaMemcpyHostToDevice); )
    
    BN254G1* outputBN254ArrayGPU;
    CUDA_CALL( cudaMalloc((void**)&outputBN254ArrayGPU, sizeof(BN254G1) * batch_size); )

    printf("launch block = %d thread = %d\n", blocks, threads_per_block);

    MSM_unit_processing <<<blocks,threads_per_block>>>( inputScalarArrayGPU, inputBaseArrayGPU, outputBN254ArrayGPU, outerc, windowSize, inner_len, batch_size);
    CUDA_CALL(cudaDeviceSynchronize());


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
 * Class:     algebra_msm_FixedBaseMSM
 * Method:    batchMSMNativeHelper
 * Signature: (IILjava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;I)[B
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_FixedBaseMSM_batchMSMNativeHelper
  (JNIEnv *env, jclass obj, jint outerc, jint windowSize, jobject multiplesOfBaseX, jobject multiplesOfBaseY, jobject multiplesOfBaseZ,  jobject bigScalars, jint BNType)
{

  jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
  jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
  jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");

  jint out_len = env->CallIntMethod(multiplesOfBaseX, java_util_ArrayList_size);
  jint inner_len = env->CallIntMethod(env->CallObjectMethod(multiplesOfBaseX, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  
  jint batch_size = env->CallIntMethod(bigScalars, java_util_ArrayList_size);

  vector<Scalar> bigScalarArray = vector<Scalar>(batch_size, Scalar());
  vector<BN254G1> multiplesOfBasePtrArray = vector<BN254G1>(out_len * inner_len, BN254G1());



  for(int i =0; i < batch_size; i++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)&bigScalarArray[i]._limbs;
      memcpy(tmp, bytes, len);
  }

  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBaseX, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray[i * inner_len + j].X._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBaseY, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray[i * inner_len + j].Y._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBaseZ, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray[i * inner_len + j].Z._limbs;
      memcpy(tmp, bytes,len);
    }
  }
  cout << "output array number of bytes=" << sizeof(BN254G1) * batch_size <<endl;
  jbyteArray resultByteArray = env->NewByteArray(sizeof(BN254G1) * batch_size);
  BN254G1* outputBN254ArrayCPU = new BN254G1[batch_size];

  fixed_batch_MSM(bigScalarArray, multiplesOfBasePtrArray, outputBN254ArrayCPU, outerc, windowSize, out_len, inner_len);
  for(int i = 0; i < batch_size; i++){
    //printMem(outputBN254ArrayCPU[i].X);
    env->SetByteArrayRegion(resultByteArray, 3 * i * sizeof(Scalar) , sizeof(Scalar) ,   reinterpret_cast<const jbyte*>(&outputBN254ArrayCPU[i].X._limbs));
    env->SetByteArrayRegion(resultByteArray, (3 * i +1) * sizeof(Scalar) , sizeof(Scalar) ,   reinterpret_cast<const jbyte*>(&outputBN254ArrayCPU[i].Y._limbs));
    env->SetByteArrayRegion(resultByteArray, (3 * i +2)* sizeof(Scalar) , sizeof(Scalar) ,   reinterpret_cast<const jbyte*>(&outputBN254ArrayCPU[i].Z._limbs));
  }

  return resultByteArray;
}




/*
 * Class:     algebra_msm_FixedBaseMSM
 * Method:    doubleBatchMSMNativeHelper
 * Signature: (IIIILjava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_FixedBaseMSM_doubleBatchMSMNativeHelper
  (JNIEnv * env, jclass obj, jint outerc1, jint windowSize1, jint outerc2, jint windowSize2, jobject multiplesOfBase1, jobject multiplesOfBase2, jobject bigScalars){

  jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
  jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
  jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");

  jint out_len1 = env->CallIntMethod(multiplesOfBase1, java_util_ArrayList_size);
  jint inner_len1 = env->CallIntMethod(env->CallObjectMethod(multiplesOfBase1, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  jint out_len2 = env->CallIntMethod(multiplesOfBase2, java_util_ArrayList_size);
  jint inner_len2 = env->CallIntMethod(env->CallObjectMethod(multiplesOfBase2, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  
  jint batch_size = env->CallIntMethod(bigScalars, java_util_ArrayList_size);
  cout << "cpp side batch size: " << batch_size << endl;


  // auto start = std::chrono::steady_clock::now();
  // vector<BigInt> bigScalarArray = vector<BigInt>(batch_size, BigInt());
  // vector<vector<BigInt>> multiplesOfBasePtrArray1 = vector<vector<BigInt>>(out_len1, vector<BigInt>(inner_len1, BigInt()));
  // vector<vector<BigInt>> multiplesOfBasePtrArray2 = vector<vector<BigInt>>(out_len2, vector<BigInt>(inner_len2, BigInt()));

  // auto end = std::chrono::steady_clock::now();
  // std::chrono::duration<double> elapsed_seconds = end-start;
  // //std::cout << "doubleBatchMSM BigInt allocation elapsed time: " << elapsed_seconds.count() << "s\n";




  // start = std::chrono::steady_clock::now();
  // for(int i =0; i < batch_size; i++){
  //     jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
  //     char* bytes = (char*)env->GetByteArrayElements(element, NULL);
  //     bigScalarArray[i].len = env->GetArrayLength(element);
  //     char* tmp = (char*)&bigScalarArray[i].bytes;

  //     memcpy(tmp +BigInt::num_of_bytes - bigScalarArray[i].len, 
  //                               bytes,
  //                               bigScalarArray[i].len);

  // }



  // //TODO parallelize it
  // for(int i = 0; i < out_len1;i++){
  //   for(int j = 0; j < inner_len1; j++){
  //     jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase1, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
  //     char* bytes = (char*)env->GetByteArrayElements(element, NULL);
  //     multiplesOfBasePtrArray1[i][j].len = env->GetArrayLength(element);
  //     char* tmp = (char*)multiplesOfBasePtrArray1[i][j].bytes;
  //     memcpy(tmp + BigInt::num_of_bytes - multiplesOfBasePtrArray1[i][j].len, bytes,  multiplesOfBasePtrArray1[i][j].len);
  //   }
  // }

  // for(int i = 0; i < out_len2;i++){
  //   for(int j = 0; j < inner_len2; j++){
  //     jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
  //     char* bytes = (char*)env->GetByteArrayElements(element, NULL);
  //     multiplesOfBasePtrArray2[i][j].len = env->GetArrayLength(element);
  //     char* tmp = (char*)multiplesOfBasePtrArray2[i][j].bytes;
  //     memcpy(tmp + BigInt::num_of_bytes - multiplesOfBasePtrArray2[i][j].len, bytes,  multiplesOfBasePtrArray2[i][j].len);    }
  // }

  //   end = std::chrono::steady_clock::now();
  //   elapsed_seconds = end-start;
  //   //std::cout << "doubleBatchMSM Read from JVM elapsed time: " << elapsed_seconds.count() << "s\n";


  // start = std::chrono::steady_clock::now();
  jbyteArray resultByteArray = env->NewByteArray( batch_size * 2);
  // for(int batch_index = 0; batch_index < batch_size; batch_index++){
  //   BigInt res1 = multiplesOfBasePtrArray1[0][0];
  //   for (int outer = 0; outer < outerc1; ++outer) {
  //       int inner = 0;
  //       for (int i = 0; i < windowSize1; ++i) {
  //           if (bigScalarArray[batch_index].testBit(outer * windowSize1 + i)) { //Returns true if and only if the designated bit is set.
  //               inner |= 1 << i;
  //           }
  //       }
  //       res1 = res1 + multiplesOfBasePtrArray1[outer][inner];
  //   }



  //   BigInt res2 = multiplesOfBasePtrArray2[0][0];
  //   for (int outer = 0; outer < outerc2; ++outer) {
  //       int inner = 0;
  //       for (int i = 0; i < windowSize2; ++i) {
  //           if (bigScalarArray[batch_index].testBit(outer * windowSize2 + i)) { //Returns true if and only if the designated bit is set.
  //               inner |= 1 << i;
  //           }
  //       }
  //       res2 = res2 + multiplesOfBasePtrArray2[outer][inner];
  //   }

  //   env->SetByteArrayRegion(resultByteArray, 2 * batch_index * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(res1.bytes));
  //   env->SetByteArrayRegion(resultByteArray, (2 * batch_index + 1) * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(res2.bytes));

  // }

  //   end = std::chrono::steady_clock::now();
  //   elapsed_seconds = end-start;
  //   //std::cout << "doubleBatchMSM C++ Compute elapsed time: " << elapsed_seconds.count() << "s\n";
  

  return resultByteArray;

  }




