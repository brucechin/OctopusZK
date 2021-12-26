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

//TODO this three Fp2 functions may have bugs
__device__ 
Fp2 add(Fp2 input1, Fp2 input2)
{
  //add should be good.
    context_t _context;
    env_t    _env(_context);  
    Fp2 result;
    memset(result.a._limbs, 0, MSM_params_t::num_of_bytes);
    memset(result.b._limbs, 0, MSM_params_t::num_of_bytes);
    Scalar modulus_binary;
    memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
    bn_t modulus;
    cgbn_load(_env, modulus, &modulus_binary);


    cgbn_add(_env, result.a, input1.a, input2.a);
    cgbn_rem(_env, result.a, result.a, modulus);

    cgbn_add(_env, result.b, input1.b, input2.b);
    cgbn_rem(_env, result.b, result.b, modulus);
    return result;
}

__device__ 
Fp2 sub(Fp2 input1, Fp2 input2)
{
    context_t _context;
    env_t    _env(_context);  
    Fp2 result;
    memset(result.a._limbs, 0, MSM_params_t::num_of_bytes);
    memset(result.b._limbs, 0, MSM_params_t::num_of_bytes);



    Scalar modulus_binary;
    memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
    bn_t modulus;
    cgbn_load(_env, modulus, &modulus_binary);



    cgbn_add(_env, result.a, input1.a, modulus);
    cgbn_sub(_env, result.a, result.a, input2.a);
    cgbn_rem(_env, result.a, result.a, modulus);

    cgbn_add(_env, result.b, input1.b, modulus);
    cgbn_sub(_env, result.b, result.b, input2.b);
    cgbn_rem(_env, result.b, result.b, modulus);
    return result;
}

__device__ 
Fp2 mul(Fp2 input1, Fp2 input2)
{
    context_t _context;
    env_t    _env(_context);  
    Fp2 result;
    memset(result.a._limbs, 0, MSM_params_t::num_of_bytes);
    memset(result.b._limbs, 0, MSM_params_t::num_of_bytes);

    Scalar modulus_binary;
    memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
    bn_t modulus;
    cgbn_load(_env, modulus, &modulus_binary);

    Scalar residue_binary;
    memcpy(residue_binary._limbs, Fp2_nonresidue_raw, MSM_params_t::num_of_bytes);
    bn_t residue;
    cgbn_load(_env, residue, &residue_binary);




    bn_t c0c0, c1c1, tmp1;

    cgbn_mul(_env, c0c0, input1.a, input2.a);
    cgbn_rem(_env, c0c0, c0c0, modulus);
    cgbn_mul(_env, c1c1, input1.b, input2.b);
    cgbn_rem(_env, c1c1, c1c1, modulus);


    cgbn_mul(_env, tmp1, residue, c1c1);
    cgbn_rem(_env, tmp1, tmp1, modulus);
    cgbn_add(_env, tmp1, tmp1, c0c0);
    cgbn_rem(_env, tmp1, tmp1, modulus);


    bn_t tmp2, tmp3;

    cgbn_add(_env, tmp2, input1.a, input1.b);
    cgbn_rem(_env, tmp2, tmp2, modulus);
    cgbn_add(_env, tmp3, input2.a, input2.b);
    cgbn_rem(_env, tmp3, tmp3, modulus);
    cgbn_mul(_env, tmp2, tmp2, tmp3);
    cgbn_rem(_env, tmp2, tmp2, modulus);
    cgbn_add(_env, tmp2, tmp2, modulus);
    cgbn_add(_env, tmp2, tmp2, modulus);
    cgbn_sub(_env, tmp2, tmp2, c0c0);
    cgbn_sub(_env, tmp2, tmp2, c1c1);
    cgbn_rem(_env, tmp2, tmp2, modulus);

    result.a = tmp1;
    result.b = tmp2;
    return result;
}


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

__device__ 
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


__device__ 
BN254G2Compute twice(BN254G2Compute a)
{
  if(isZero(a)){
    return a;
  }

  BN254G2Compute result;
  memset(result.X.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.X.b._limbs, 0, MSM_params_t::num_of_bytes);

  memset(result.Y.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Y.b._limbs, 0, MSM_params_t::num_of_bytes);

  memset(result.Z.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Z.b._limbs, 0, MSM_params_t::num_of_bytes);

  Fp2 a_x, a_y, a_z;
  a_x = a.X;
  a_y = a.Y;
  a_z = a.Z;


  Fp2 A,B,C,D,E,F,X3,Y3,Z3, eightC;


  A = mul(a_x, a_x);
  B = mul(a_y, a_y);

  C = mul(B, B);

 // D = 2 * ((X1 + B)^2 - A - C)
  D = add(a_x, B);
  D = mul(D, D);
  D = sub(D, A);
  D = sub(D, C);
  D = add(D, D);

  // E = 3 * A
  E = add(A, A);
  E = add(E, A);

  // F = E^2
  F = mul(E, E);

   // X3 = F - 2 D
  X3 = sub(F, D);
  X3 = sub(X3, D);



  eightC = add(C, C);
  eightC = add(eightC, eightC);
  eightC = add(eightC, eightC);

  // Y3 = E * (D - X3) - 8 * C
  Y3 = sub(D, X3);
  Y3 = mul(E, Y3);
  Y3 = sub(Y3, eightC);


  // Z3 = 2 * Y1 * Z1
  Z3 = mul(a_y, a_z);
  Z3 = add(Z3, Z3);

  result.X = X3;
  result.Y = Y3;
  result.Z = Z3;

  return result;

}


void printMem(Scalar input)
{
    for(int i = 0; i < MSM_params_t::BITS/32; i++){
      std::bitset<32> tmp(input._limbs[i]);
      std::cout << tmp << "|";
    }
    printf("finished\n");
}


__device__ 
void print_bn_t(bn_t &number, int instance_id_) {
  using __env_t = bn_t::parent_env_t;
  const int IPB = 128/__env_t::TPI;
  const int TPI = __env_t::TPI;
  __shared__ uint32_t n[IPB][(__env_t::BITS/32)] ;
  __shared__ uint32_t vote[IPB];
  bool is_represent = (threadIdx.x % TPI) == 0;
  int  instance_id  = threadIdx.x / TPI;
  int  tid_in_instance = threadIdx.x % TPI;
  int global_instance_id = (threadIdx.x + blockIdx.x * blockDim.x)/TPI ;
if(global_instance_id == instance_id_){
  if (is_represent) vote[instance_id] = 0;
  for (int i = 0; i < __env_t::LIMBS; i++)
    n[instance_id][tid_in_instance * __env_t::LIMBS + i] = number._limbs[i];
  atomicAdd(&vote[instance_id], 1);
  while (vote[instance_id] < TPI) ;
  if (is_represent) {
    //printf("instance %d is ", global_instance_id);
    for (int i = 0; i < __env_t::BITS/32; i++) {
      printf(" %u |", n[instance_id][i]);
    }
    printf("\n");
  }
}

  
}

__device__ 
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

__device__ 
BN254G2Compute add(BN254G2Compute a, BN254G2Compute b) {
  // Handle special cases having to do with O

  if (isZero(a)) {
      return b;
  }

  if (isZero(b)) {
      return a;
  }

  context_t _context;
  env_t    _env(_context);

  BN254G2Compute result;
  memset(result.X.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.X.b._limbs, 0, MSM_params_t::num_of_bytes);

  memset(result.Y.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Y.b._limbs, 0, MSM_params_t::num_of_bytes);

  memset(result.Z.a._limbs, 0, MSM_params_t::num_of_bytes);
  memset(result.Z.b._limbs, 0, MSM_params_t::num_of_bytes);

  Fp2 a_x, a_y, a_z, b_x, b_y, b_z;
  a_x = a.X; 
  a_y = a.Y; 
  a_z = a.Z;
  b_x = b.X;
  b_y = b.Y;
  b_z = b.Z;

  Fp2 Z1Z1, Z2Z2, U1, U2, Z1_cubed, Z2_cubed, S1, S2;
  
  Z1Z1 = mul(a_z, a_z);
  Z2Z2 = mul(b_z, b_z);

  U1 = mul(a_x, Z2Z2);
  U2 = mul(b_x, Z1Z1);

  Z1_cubed = mul(a_z, Z1Z1);
  Z2_cubed = mul(b_z, Z2Z2);

  S1 = mul(a_y, Z2_cubed);
  S2 = mul(b_y, Z1_cubed);


  if (cgbn_equals(_env, U1.a, U2.a) && cgbn_equals(_env, U1.b, U2.b)
    && cgbn_equals(_env, S1.a, S2.a) && cgbn_equals(_env, S1.b, S2.b)
    ) {
      // Double case; nothing above can be reused.
      return twice(a);
  }

  Fp2 H, S2_minus_S1, I, J, r, V, X3, S1_J, Y3, Z3;
  
  // H = U2-U1
  H = sub(U2, U1);
  S2_minus_S1 = sub(S2, S1);

  // I = (2 * H)^2
  I = add(H, H);
  I = mul(I, I);

  // J = H * I
  J = mul(H, I);

  // r = 2 * (S2-S1)
  r = add(S2_minus_S1, S2_minus_S1);

  // V = U1 * I
  V = mul(U1, I);

  // X3 = r^2 - J - 2 * V
  X3 = mul(r, r);
  X3 = sub(X3, J);
  X3 = sub(X3, V);
  X3 = sub(X3, V);

  // Y3 = r * (V-X3)-2 * S1_J
  S1_J = mul(S1, J);
  Y3 = sub(V, X3);
  Y3 = mul(r, Y3);
  Y3 = sub(Y3, S1_J);
  Y3 = sub(Y3, S1_J);

  // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

  Z3 = add(a_z, b_z);
  Z3 = mul(Z3, Z3);
  Z3 = sub(Z3, Z1Z1);
  Z3 = sub(Z3, Z2Z2);
  Z3 = mul(Z3, H);

  result.X = X3;
  result.Y = Y3;
  result.Z = Z3;

  return result;

}


__global__ void fixedbase_MSM_unit_processing_G1(Scalar* inputScalarArray, BN254G1* inputBaseArray, BN254G1* outputBN254Array, int outerc, int windowSize, int tableInnerSize, int batch_size){
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
    for (int outer = 0; outer < outerc; ++outer) {
        int inner = 0;
        for (int i = 0; i < windowSize; ++i) {
            //testBit is correct
            if (testBit(inputScalarArray[idx], outer * windowSize + i)) {
                inner |= 1 << i;
            }
        }
        //printf("idx=%d, outer=%d, inner=%d\n", idx, outer, inner);
        BN254G1Compute to_added;
        cgbn_load(_env, to_added.X, &inputBaseArray[tableInnerSize * outer + inner].X);
        cgbn_load(_env, to_added.Y, &inputBaseArray[tableInnerSize * outer + inner].Y);
        cgbn_load(_env, to_added.Z, &inputBaseArray[tableInnerSize * outer + inner].Z);
        res = add(res, to_added);
    }

    cgbn_store(_env, &outputBN254Array[idx].X, res.X);
    cgbn_store(_env, &outputBN254Array[idx].Y, res.Y);
    cgbn_store(_env, &outputBN254Array[idx].Z, res.Z);

    return;
}


__global__ void fixedbase_MSM_unit_processing_G2(Scalar* inputScalarArray, BN254G2* inputBaseArray, BN254G2* outputBN254Array, int outerc, int windowSize, int tableInnerSize, int batch_size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    context_t _context;
    env_t    _env(_context);
    BN254G2Compute res;
    cgbn_load(_env, res.X.a, &inputBaseArray[0].Xa);
    cgbn_load(_env, res.X.b, &inputBaseArray[0].Xb);
    cgbn_load(_env, res.Y.a, &inputBaseArray[0].Ya);
    cgbn_load(_env, res.Y.b, &inputBaseArray[0].Yb);
    cgbn_load(_env, res.Z.a, &inputBaseArray[0].Za);
    cgbn_load(_env, res.Z.b, &inputBaseArray[0].Zb);


    if(idx >= batch_size){
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
        //printf("idx=%d, outer=%d, inner=%d\n", idx, outer, inner);
        BN254G2Compute to_added;
        cgbn_load(_env, to_added.X.a, &inputBaseArray[tableInnerSize * outer + inner].Xa);
        cgbn_load(_env, to_added.X.b, &inputBaseArray[tableInnerSize * outer + inner].Xb);
        cgbn_load(_env, to_added.Y.a, &inputBaseArray[tableInnerSize * outer + inner].Ya);
        cgbn_load(_env, to_added.Y.b, &inputBaseArray[tableInnerSize * outer + inner].Yb);
        cgbn_load(_env, to_added.Z.a, &inputBaseArray[tableInnerSize * outer + inner].Za);
        cgbn_load(_env, to_added.Z.b, &inputBaseArray[tableInnerSize * outer + inner].Zb);

        res = add(res, to_added);
    }

    cgbn_store(_env, &outputBN254Array[idx].Xa, res.X.a);
    cgbn_store(_env, &outputBN254Array[idx].Xb, res.X.b);
    cgbn_store(_env, &outputBN254Array[idx].Ya, res.Y.a);
    cgbn_store(_env, &outputBN254Array[idx].Yb, res.Y.b);
    cgbn_store(_env, &outputBN254Array[idx].Za, res.Z.a);
    cgbn_store(_env, &outputBN254Array[idx].Zb, res.Z.b);


    return;
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
    CUDA_CALL( cudaMemset(outputBN254ArrayGPU, 0, sizeof(BN254G1) * batch_size); )

    printf("launch block = %d thread = %d\n", blocks, threads_per_block);

    fixedbase_MSM_unit_processing_G1 <<<blocks,threads_per_block, 32 * 1024>>>( inputScalarArrayGPU, inputBaseArrayGPU, outputBN254ArrayGPU, outerc, windowSize, inner_len, batch_size);
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


void  fixed_double_batch_MSM(std::vector<Scalar> & bigScalarArray, std::vector<BN254G1> &multiplesOfBase1PtrArray, std::vector<BN254G2> &multiplesOfBase2PtrArray,
                            BN254G1* outputArrayG1, BN254G2* outputArrayG2,
                            int outerc1, int windowSize1, int outerc2, int windowSize2, 
                            int out_len1, int inner_len1 , int out_len2, int inner_len2)
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
    
    BN254G1 *inputBase1ArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputBase1ArrayGPU, sizeof(BN254G1) * out_len1 * inner_len1); )
    CUDA_CALL( cudaMemcpy(inputBase1ArrayGPU, (void**)&multiplesOfBase1PtrArray[0], sizeof(BN254G1) * out_len1 * inner_len1, cudaMemcpyHostToDevice); )
    
    BN254G2 *inputBase2ArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputBase2ArrayGPU, sizeof(BN254G2) * out_len2 * inner_len2); )
    CUDA_CALL( cudaMemcpy(inputBase2ArrayGPU, (void**)&multiplesOfBase2PtrArray[0], sizeof(BN254G2) * out_len2 * inner_len2, cudaMemcpyHostToDevice); )
    

    BN254G1* outputBN254G1ArrayGPU;
    CUDA_CALL( cudaMalloc((void**)&outputBN254G1ArrayGPU, sizeof(BN254G1) * batch_size); )
    CUDA_CALL( cudaMemset(outputBN254G1ArrayGPU, 0, sizeof(BN254G1) * batch_size); )

    BN254G2* outputBN254G2ArrayGPU;
    CUDA_CALL( cudaMalloc((void**)&outputBN254G2ArrayGPU, sizeof(BN254G2) * batch_size); )
    CUDA_CALL( cudaMemset(outputBN254G2ArrayGPU, 0, sizeof(BN254G2) * batch_size); )

    printf("launch block = %d thread = %d\n", blocks, threads_per_block);
    fixedbase_MSM_unit_processing_G1 <<<blocks,threads_per_block, 32 * 1024>>>( inputScalarArrayGPU, inputBase1ArrayGPU, outputBN254G1ArrayGPU, outerc1, windowSize1, inner_len1, batch_size);
    
    
    CUDA_CALL(cudaDeviceSynchronize());

    fixedbase_MSM_unit_processing_G2 <<<blocks,threads_per_block, 32 * 1024>>> (inputScalarArrayGPU, inputBase2ArrayGPU, outputBN254G2ArrayGPU, outerc2, windowSize2, inner_len2, batch_size);

    CUDA_CALL(cudaDeviceSynchronize());


    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }

    CUDA_CALL(cudaMemcpy((void**)outputArrayG1, outputBN254G1ArrayGPU, sizeof(BN254G1) * batch_size, cudaMemcpyDeviceToHost); )
    CUDA_CALL(cudaMemcpy((void**)outputArrayG2, outputBN254G2ArrayGPU, sizeof(BN254G2) * batch_size, cudaMemcpyDeviceToHost); )

    CUDA_CALL(cudaDeviceSynchronize());

    CUDA_CALL(cudaFree(inputScalarArrayGPU));
    CUDA_CALL(cudaFree(inputBase1ArrayGPU));
    CUDA_CALL(cudaFree(outputBN254G1ArrayGPU));
    CUDA_CALL(cudaFree(inputBase2ArrayGPU));
    CUDA_CALL(cudaFree(outputBN254G2ArrayGPU));
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

  cout << "CUDA side base out and inner len :" << out_len << " " << inner_len <<endl;

  cout << "CUDA side base outerc and windowSize :" << outerc << " " << windowSize <<endl;

  auto start = std::chrono::steady_clock::now();

  for(int i =0; i < batch_size; i++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      //cout << len << " ";
      char* tmp = (char*)&bigScalarArray[i]._limbs;
     memcpy(tmp, bytes, len);
  }

  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBaseX, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
           // cout << len << " ";

      char* tmp = (char*)multiplesOfBasePtrArray[i * inner_len + j].X._limbs;
     memcpy(tmp, bytes,len);

    }
  }

  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBaseY, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      //cout << len << " ";

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
  memset(outputBN254ArrayCPU, 0, sizeof(BN254G1) * batch_size);

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "C++ read from JVM elapsed time: " << elapsed_seconds.count() << "s\n";

  fixed_batch_MSM(bigScalarArray, multiplesOfBasePtrArray, outputBN254ArrayCPU, outerc, windowSize, out_len, inner_len);
  end = std::chrono::steady_clock::now();


  start = std::chrono::steady_clock::now();
  for(int i = 0; i < batch_size; i++){
    env->SetByteArrayRegion(resultByteArray, i * sizeof(BN254G1) , sizeof(BN254G1) ,   reinterpret_cast<const jbyte*>(&outputBN254ArrayCPU[i]));
  }
  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "C++ write results back to JVM elapsed time: " << elapsed_seconds.count() << "s\n";


  return resultByteArray;
}




/*
 * Class:     algebra_msm_FixedBaseMSM
 * Method:    doubleBatchMSMNativeHelper
 * Signature: (IIIILjava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */

JNIEXPORT jbyteArray JNICALL Java_algebra_msm_FixedBaseMSM_doubleBatchMSMNativeHelper
  (JNIEnv * env, jclass obj, jint outerc1, jint windowSize1, jint outerc2, jint windowSize2, jobject multiplesOfBase1_X, jobject multiplesOfBase1_Y,
   jobject multiplesOfBase1_Z, jobject multiplesOfBase2_Xa, jobject multiplesOfBase2_Ya, jobject multiplesOfBase2_Za, jobject multiplesOfBase2_Xb,
    jobject multiplesOfBase2_Yb, jobject multiplesOfBase2_Zb, jobject bigScalars){

  jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
  jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
  jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");

  jint out_len1 = env->CallIntMethod(multiplesOfBase1_X, java_util_ArrayList_size);
  jint inner_len1 = env->CallIntMethod(env->CallObjectMethod(multiplesOfBase1_X, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  jint out_len2 = env->CallIntMethod(multiplesOfBase2_Xa, java_util_ArrayList_size);
  jint inner_len2 = env->CallIntMethod(env->CallObjectMethod(multiplesOfBase2_Xa, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  
  jint batch_size = env->CallIntMethod(bigScalars, java_util_ArrayList_size);
  cout << "cpp side batch size: " << batch_size << endl;
  cout << "cpp side out len1 , inner_len1 : " << out_len1 << " " << inner_len1 << endl;
  cout << "cpp side out len2 , inner_len2 : " << out_len2 << " " << inner_len2 << endl;

  auto start = std::chrono::steady_clock::now();
  vector<Scalar> bigScalarArray = vector<Scalar>(batch_size, Scalar());
  vector<BN254G1> multiplesOfBase1PtrArray = vector<BN254G1>(out_len1 * inner_len1, BN254G1());
  vector<BN254G2> multiplesOfBase2PtrArray = vector<BN254G2>(out_len2 * inner_len2, BN254G2());

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "doubleBatchMSM BigInt allocation elapsed time: " << elapsed_seconds.count() << "s\n";




  start = std::chrono::steady_clock::now();
  for(int i =0; i < batch_size; i++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)&bigScalarArray[i]._limbs;
      memcpy(tmp, bytes, len);
  }
  //std::cout << "biscalar" <<endl;

  for(int i = 0; i < out_len1;i++){
    for(int j = 0; j < inner_len1; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase1_X, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase1PtrArray[i * inner_len1 + j].X._limbs;
      memcpy(tmp, bytes,len);

    }
  }

  for(int i = 0; i < out_len1;i++){
    for(int j = 0; j < inner_len1; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase1_Y, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase1PtrArray[i * inner_len1 + j].Y._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len1;i++){
    for(int j = 0; j < inner_len1; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase1_Z, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase1PtrArray[i * inner_len1 + j].Z._limbs;
      memcpy(tmp, bytes,len);
    }
  }
  //std::cout << "G1" <<endl;



  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Xa, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Xa._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Xb, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Xb._limbs;
      memcpy(tmp, bytes,len);
    }
  }
  //std::cout << "G2X" <<endl;

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Ya, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Ya._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Yb, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Yb._limbs;
      memcpy(tmp, bytes,len);
    }
  }
  //std::cout << "G2Y" <<endl;

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Za, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Za._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2_Zb, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      int len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBase2PtrArray[i * inner_len2 + j].Zb._limbs;
      memcpy(tmp, bytes,len);
    }
  }

  //std::cout << "G2" <<endl;

  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "doubleBatchMSM Read from JVM elapsed time: " << elapsed_seconds.count() << "s\n";


  jbyteArray resultByteArray = env->NewByteArray(batch_size * (sizeof(BN254G1) + sizeof(BN254G2)));
  BN254G1* outputBN254G1ArrayCPU = new BN254G1[batch_size];
  memset(outputBN254G1ArrayCPU, 0, sizeof(BN254G1) * batch_size);
  BN254G2* outputBN254G2ArrayCPU = new BN254G2[batch_size];
  memset(outputBN254G2ArrayCPU, 0, sizeof(BN254G2) * batch_size);

  fixed_double_batch_MSM(bigScalarArray, multiplesOfBase1PtrArray, multiplesOfBase2PtrArray, 
                          outputBN254G1ArrayCPU, outputBN254G2ArrayCPU,
                          outerc1, windowSize1, outerc2, windowSize2,
                          out_len1, inner_len1, out_len2, inner_len2);

  start = std::chrono::steady_clock::now();

  for(int i = 0; i < batch_size; i++){
    env->SetByteArrayRegion(resultByteArray, i * (sizeof(BN254G1) +sizeof(BN254G2)), sizeof(BN254G1) ,   reinterpret_cast<const jbyte*>(&outputBN254G1ArrayCPU[i]));
    env->SetByteArrayRegion(resultByteArray, i * (sizeof(BN254G1) +sizeof(BN254G2)) + sizeof(BN254G1), sizeof(BN254G2) ,   reinterpret_cast<const jbyte*>(&outputBN254G2ArrayCPU[i]));
  }

  end = std::chrono::steady_clock::now();
  elapsed_seconds = end-start;
  std::cout << "doubleBatchMSM C++ Compute elapsed time: " << elapsed_seconds.count() << "s\n";
  
  return resultByteArray;

  }




















// __device__ __forceinline__
// bool equals(BN254G1 a, BN254G1 b)
// {
//   if(isZero(a)){
//     return isZero(b);
//   }

//   if(isZero(b)){
//     return false;
//   }

//   context_t _context;
//   env_t    _env(_context);

//   Scalar modulus_binary;
//   memcpy(modulus_binary._limbs, modulus_raw_G1, MSM_params_t::num_of_bytes);
//   bn_t modulus;
//   cgbn_load(_env, modulus, &modulus_binary);

//   bn_t a_x, a_y, a_z, b_x, b_y, b_z;
//   cgbn_load(_env, a_x, &a.X);
//   cgbn_load(_env, a_y, &a.Y);
//   cgbn_load(_env, a_z, &a.Z);
//   cgbn_load(_env, b_x, &b.X);
//   cgbn_load(_env, b_y, &b.Y);
//   cgbn_load(_env, b_z, &b.Z);   


//   bn_t Z1_squared, Z2_squared, XZ1_squared, XZ2_squared, Z1_cubed, Z2_cubed;
//   cgbn_mul(_env, Z1_squared, a_z, a_z);
//   cgbn_rem(_env, Z1_squared, Z1_squared, modulus);

//   cgbn_mul(_env, Z2_squared, b_z, b_z);
//   cgbn_rem(_env, Z2_squared, Z2_squared, modulus);

//   cgbn_mul(_env, XZ1_squared, Z1_squared, b_x);
//   cgbn_rem(_env, XZ1_squared, XZ1_squared, modulus);

//   cgbn_mul(_env, XZ2_squared, Z2_squared, a_x);
//   cgbn_rem(_env, XZ2_squared, XZ2_squared, modulus);


//   if(cgbn_equals(_env, XZ1_squared, XZ2_squared)){
//     return false;
//   }

//   cgbn_mul(_env, Z1_cubed, Z1_squared, a_z);
//   cgbn_rem(_env, Z1_cubed, Z1_cubed, modulus);

//   cgbn_mul(_env, Z2_cubed, Z2_squared, b_z);
//   cgbn_rem(_env, Z2_cubed, Z2_cubed, modulus);


//   bn_t YZ2_cubed, YZ1_cubed;
//   cgbn_mul(_env, YZ1_cubed, Z1_cubed, b_y);
//   cgbn_rem(_env, YZ1_cubed, YZ1_cubed, modulus);

//   cgbn_mul(_env, YZ2_cubed, Z2_cubed, a_y);
//   cgbn_rem(_env, YZ2_cubed, YZ2_cubed, modulus);

//   if(cgbn_equals(_env, YZ1_cubed, YZ2_cubed)){
//     return false;
//   }

//   return true;
// }




// this code is only used for Fp2 debugging purpose.
// __global__ void fp2_test(BN254G2* inputBaseArray, BN254G2* outputArray){
//       context_t _context;
//     env_t    _env(_context);
//         const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
//     Fp2 a, b;


//     cgbn_load(_env, a.a, &inputBaseArray[idx].Xa);
//     cgbn_load(_env, a.b, &inputBaseArray[idx].Xb);
//     cgbn_load(_env, b.a, &inputBaseArray[idx].Ya);
//     cgbn_load(_env, b.b, &inputBaseArray[idx].Yb);
//     // print_bn_t(a.a, 1);
//     // print_bn_t(a.b, 1);
//     // print_bn_t(b.a, 1);
//     // print_bn_t(b.b, 1);
//     Fp2 a_b_add = add(a, b);
//     Fp2 a_b_sub = sub(a, b);
//     Fp2 a_b_mul = mul(a, b);

//     // print_bn_t(a_b_add.a, 1);
//     // print_bn_t(a_b_add.b, 1);
//     cgbn_store(_env, &outputArray[idx].Xa, a_b_add.a);
//     cgbn_store(_env, &outputArray[idx].Xb, a_b_add.b);
//     cgbn_store(_env, &outputArray[idx].Ya, a_b_sub.a);
//     cgbn_store(_env, &outputArray[idx].Yb, a_b_sub.b);
//     cgbn_store(_env, &outputArray[idx].Za, a_b_mul.a);
//     cgbn_store(_env, &outputArray[idx].Zb, a_b_mul.b);

// }