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

void printMem(Scalar input)
{
    for(int i = 0; i < MSM_params_t::BITS/32; i++){
      std::cout << input._limbs[i] << "|";
    }
    printf("finished\n");
}

__global__ void vector_print(int* input, int size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    if(idx == 0){
        for(int i = 0; i < size; i++) {
            printf("%d|", input[i]);
        }
    }
}

__global__ void prefix_sum_first_step(int* input, int size, int stride){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    int index = (idx + 1) * stride * 2 -1;
    if(index < size){
        input[index] += input[index - stride];
    }
}

__global__ void prefix_sum_second_step(int* input, int size, int stride){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    int index = (idx + 1) *stride *2 -1;
    if(index + stride <size){
        input[index + stride] += input[index];
    }
}


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


__device__ 
bool testBit(Scalar input, int n)
{
    int byte_index =n / 32;
    int byte_offset = n % 32;
    return CHECK_BIT(input._limbs[byte_index], byte_offset);
}

__device__ 
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

__device__ 
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
  cgbn_rem(_env, D, D, modulus);
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
  cgbn_add(_env, X3, X3, modulus);
  cgbn_sub(_env, X3, X3, D);
  cgbn_sub(_env, X3, X3, D);
  cgbn_rem(_env, X3, X3, modulus);




  cgbn_add(_env, eightC, C, C);
  cgbn_rem(_env, eightC, eightC, modulus);
  cgbn_add(_env, eightC, eightC, eightC);
  cgbn_rem(_env, eightC, eightC, modulus);
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


  //print_bn_t(X3, 1);
  result.X = X3;
  result.Y = Y3;
  result.Z = Z3;

  return result;

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
      //printf("twice is called");
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



__device__ void printMemGPU(Scalar input)
{
    for(int i = 0; i < MSM_params_t::BITS/32; i++){
      printf("%u|", input._limbs[i] );
    }
    printf("finished\n");
}


__global__ void pippengerMSM_unit1(Scalar* inputScalarArray,  int* inputBucketMappingLocation, int* bucketCounter, int* bucketIndex,  int batch_size, int c, int k, int numGroups){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    if(idx >= batch_size){
        return;
    }

    int id = 0;
    for (int j = 0; j < c; j++) {
        if (testBit(inputScalarArray[idx], k * c + j)) {
            id |= 1 << j;
        }
    }

    int one = 1;
    bucketIndex[idx] = atomicAdd((bucketCounter + id), one);
    inputBucketMappingLocation[idx] = id;

    return;
}

__global__ void pippengerMSMG1_unit2_write_together(BN254G1* inputBaseArray, BN254G1* outputBaseArray, int* inputBucketMappingLocation,  int* bucketCounter, int* bucketIndex, int batch_size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);
    context_t _context;
    env_t    _env(_context);
    
    if(idx >= batch_size){
      return;
    }

    int outputIndex;
    if(inputBucketMappingLocation[idx] == 0){
        outputIndex = bucketIndex[idx];
    }else{
        outputIndex = bucketCounter[inputBucketMappingLocation[idx] - 1] + bucketIndex[idx];
    }

    outputBaseArray[outputIndex] = inputBaseArray[idx];
    return;
}

__global__ void pippengerMSMG1_unit2_reduce_to_buckets(BN254G1* outputBaseArray, BN254G1* outputBuckets,  int* bucketCounter, int numBuckets, BN254G1* zero){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    context_t _context;
    env_t    _env(_context);
    
    if(idx >= numBuckets){
      return;
    }

    int left, right; 
    if(idx == 0){
        return ; //when bucketID = 0, skip.
    }else{
        left = bucketCounter[idx - 1];
        right = bucketCounter[idx];
    }

    BN254G1Compute res;
    cgbn_load(_env, res.X, &zero[0].X);
    cgbn_load(_env, res.Y, &zero[0].Y);
    cgbn_load(_env, res.Z, &zero[0].Z);

    for(; left < right; left++){
        BN254G1Compute to_add; 
        cgbn_load(_env, to_add.X, &outputBaseArray[left].X);
        cgbn_load(_env, to_add.Y, &outputBaseArray[left].Y);
        cgbn_load(_env, to_add.Z, &outputBaseArray[left].Z);

        res = add(res, to_add);
    }

    cgbn_store(_env, &outputBuckets[idx].X, res.X);
    cgbn_store(_env, &outputBuckets[idx].Y, res.Y);
    cgbn_store(_env, &outputBuckets[idx].Z, res.Z);
    return;
}

__global__ void prefix_sum_G1_reverse_first_step(BN254G1* buckets, int size, int stride){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    int index = idx * stride * 2 ;
    context_t _context;
    env_t    _env(_context);

    if(index + stride < size){
        BN254G1Compute res, to_add;
        cgbn_load(_env, res.X, &buckets[index].X);
        cgbn_load(_env, res.Y, &buckets[index].Y);
        cgbn_load(_env, res.Z, &buckets[index].Z);
        cgbn_load(_env, to_add.X, &buckets[index + stride].X);
        cgbn_load(_env, to_add.Y, &buckets[index + stride].Y);
        cgbn_load(_env, to_add.Z, &buckets[index + stride].Z);
        res = add(res, to_add);
        cgbn_store(_env, &buckets[index].X, res.X);
        cgbn_store(_env, &buckets[index].Y, res.Y);
        cgbn_store(_env, &buckets[index].Z, res.Z);
    }

}

__global__ void prefix_sum_G1_reverse_second_step(BN254G1* buckets, BN254G1* result, int size, int stride, bool copyToResult){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    int index = idx *stride *2 ;
    context_t _context;
    env_t    _env(_context);
    if(index <size &&index > 0){
        BN254G1Compute res, to_add;
        cgbn_load(_env, res.X, &buckets[index - stride].X);
        cgbn_load(_env, res.Y, &buckets[index - stride].Y);
        cgbn_load(_env, res.Z, &buckets[index - stride].Z);
        cgbn_load(_env, to_add.X, &buckets[index ].X);
        cgbn_load(_env, to_add.Y, &buckets[index ].Y);
        cgbn_load(_env, to_add.Z, &buckets[index ].Z);
        res = add(res, to_add);
        cgbn_store(_env, &buckets[index  - stride].X, res.X);
        cgbn_store(_env, &buckets[index  - stride].Y, res.Y);
        cgbn_store(_env, &buckets[index  - stride].Z, res.Z);

    }

    
}

__global__ void pippengerMSMG1_unit2_final_add(BN254G1* result, BN254G1* buckets){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    context_t _context;
    env_t    _env(_context);
    if(idx == 0){
                BN254G1Compute res,to_add;
        cgbn_load(_env, res.X, &result[0].X);
        cgbn_load(_env, res.Y, &result[0].Y);
        cgbn_load(_env, res.Z, &result[0].Z);
        cgbn_load(_env, to_add.X, &buckets[1].X);
        cgbn_load(_env, to_add.Y, &buckets[1].Y);
        cgbn_load(_env, to_add.Z, &buckets[1].Z);
        res = add(res, to_add);
        cgbn_store(_env, &result[0].X, res.X);
        cgbn_store(_env, &result[0].Y, res.Y);
        cgbn_store(_env, &result[0].Z, res.Z); 
    }
}

__global__ void pippengerMSMG1_unit3(BN254G1* result, int c){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x)/MSM_params_t::TPI;
    context_t _context;
    env_t    _env(_context);


    //because c is usually very small, we do not need to parallelize this step.
    if(idx == 0){
        BN254G1Compute res;
        cgbn_load(_env, res.X, &result[0].X);
        cgbn_load(_env, res.Y, &result[0].Y);
        cgbn_load(_env, res.Z, &result[0].Z);
        for (int i = 0; i < c; i++) {
            res = twice(res);
        } 
        cgbn_store(_env, &result[0].X, res.X);
        cgbn_store(_env, &result[0].Y, res.Y);
        cgbn_store(_env, &result[0].Z, res.Z);

    }


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
    CUDA_CALL(cudaSetDevice(0));
    printf("finish CUDA device set\n");

    //TODO lianke: sometimes get stuck in memcpy.
    Scalar *inputScalarArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputScalarArrayGPU, sizeof(Scalar) * batch_size); )
    CUDA_CALL( cudaMemcpy(inputScalarArrayGPU, (void**)&bigScalarArray[0], sizeof(Scalar) * batch_size, cudaMemcpyHostToDevice); )
    
    BN254G1 *inputBaseArrayGPU; 
    CUDA_CALL( cudaMalloc((void**)&inputBaseArrayGPU, sizeof(BN254G1) * batch_size); )
    CUDA_CALL( cudaMemcpy(inputBaseArrayGPU, (void**)&multiplesOfBasePtrArray[0], sizeof(BN254G1) * batch_size, cudaMemcpyHostToDevice); )
    CUDA_CALL(cudaDeviceSynchronize());

    int numBits = 254;//BN254 specific value;
    int length = batch_size;
    int log2Length = max(1, (int)log2(length));
    int c= log2Length - (log2Length/3);
    int numBuckets = 1 << c;
    int numGroups = (numBits + c - 1 ) / c;
    BN254G1  zero;
    memset(&zero, 0, sizeof(BN254G1));
    zero.Y._limbs[0] = 1;
    vector<BN254G1> buketsModel(numBuckets, zero);
    BN254G1* resultGPU ;
    CUDA_CALL( cudaMalloc((void**)&resultGPU, sizeof(BN254G1)); )
    CUDA_CALL( cudaMemcpy(resultGPU, &zero, sizeof(BN254G1), cudaMemcpyHostToDevice); )
    BN254G1* zeroGPU ;
    CUDA_CALL( cudaMalloc((void**)&zeroGPU, sizeof(BN254G1)); )
    CUDA_CALL( cudaMemcpy(zeroGPU, &zero, sizeof(BN254G1), cudaMemcpyHostToDevice); )
    CUDA_CALL(cudaDeviceSynchronize());

    for(int k = numGroups - 1; k >= 0; k--){
        printf("k=%d ", k);
        int* bucketCounter;
        int* bucketIndex;
        int* inputBucketMappingLocation;
        BN254G1* bucketsBeforeAggregate;
        CUDA_CALL( cudaMalloc((void**)&bucketCounter, sizeof(int) * numBuckets);)
        CUDA_CALL( cudaMemset(bucketCounter, 0, numBuckets * sizeof(int));)
        CUDA_CALL( cudaMalloc((void**)&bucketIndex, sizeof(int) * batch_size);)
        CUDA_CALL( cudaMemset(bucketIndex, 0, batch_size  * sizeof(int));)
        CUDA_CALL( cudaMalloc((void**)&inputBucketMappingLocation, sizeof(int) * batch_size);)
        CUDA_CALL( cudaMemset(inputBucketMappingLocation, 0, batch_size  * sizeof(int));)
        CUDA_CALL(cudaMalloc((void**)&bucketsBeforeAggregate, sizeof(BN254G1) * batch_size);)
        CUDA_CALL( cudaMemset(bucketsBeforeAggregate, 0, sizeof(BN254G1) * batch_size);)

        BN254G1* buckets;
        CUDA_CALL( cudaMalloc((void**)&buckets, sizeof(BN254G1) * numBuckets); )
        CUDA_CALL(cudaMemcpy(buckets, (void**)&buketsModel[0], sizeof(BN254G1) * numBuckets, cudaMemcpyHostToDevice);)
        int num_blocks_unit1 = (batch_size + threads_per_block - 1) / threads_per_block;
        pippengerMSM_unit1 <<<num_blocks_unit1,threads_per_block>>>( inputScalarArrayGPU, inputBucketMappingLocation, bucketCounter, bucketIndex,  batch_size, c, k, numGroups);
        CUDA_CALL(cudaDeviceSynchronize();)


        int num_blocks_prefix_sum = (numBuckets + threads_per_block - 1)/threads_per_block;
        for(int stride = 1; stride <= numBuckets/2; stride = stride *2){
            prefix_sum_first_step<<< num_blocks_prefix_sum, threads_per_block>>>(bucketCounter, numBuckets, stride);
            CUDA_CALL(cudaDeviceSynchronize();)
        }
        
        for(int stride = numBuckets/4; stride > 0; stride /= 2){
            prefix_sum_second_step<<< num_blocks_prefix_sum, threads_per_block>>>(bucketCounter, numBuckets, stride);
            CUDA_CALL(cudaDeviceSynchronize();)
        }
        

        //now bucketCounter is prefix-sumed version. write into bucketsBeforeAggregate
        int num_blocks_write_together = (batch_size + threads_per_block - 1)/threads_per_block;
        pippengerMSMG1_unit2_write_together<<<num_blocks_write_together, threads_per_block>>>(inputBaseArrayGPU, bucketsBeforeAggregate, inputBucketMappingLocation, bucketCounter, bucketIndex, batch_size);
        CUDA_CALL(cudaDeviceSynchronize();)

        //reduce the bucketsBeforeAggregate, because they have been write to ajacent positions.
        int num_blocks_reduce_buckets = (numBuckets + instance_per_block - 1)/instance_per_block;
        pippengerMSMG1_unit2_reduce_to_buckets<<<num_blocks_reduce_buckets, threads_per_block>>>(bucketsBeforeAggregate, buckets, bucketCounter, numBuckets, zeroGPU);
        CUDA_CALL(cudaDeviceSynchronize();)


        int num_blocks_prefix_sum_G1 = (numBuckets + instance_per_block - 1)/instance_per_block;


        for(int stride = 1; stride <= numBuckets/2; stride = stride *2){
            prefix_sum_G1_reverse_first_step<<< num_blocks_prefix_sum_G1, threads_per_block>>>(buckets, numBuckets, stride);
            CUDA_CALL(cudaDeviceSynchronize();)
        }
        
        for(int stride = numBuckets/4; stride > 0; stride /= 2){
            prefix_sum_G1_reverse_second_step<<< num_blocks_prefix_sum_G1, threads_per_block>>>(buckets, resultGPU, numBuckets, stride, false);
            CUDA_CALL(cudaDeviceSynchronize();)
        }
        
        for(int stride = 1; stride <= numBuckets/2; stride = stride *2){
            prefix_sum_G1_reverse_first_step<<< num_blocks_prefix_sum_G1, threads_per_block>>>(buckets, numBuckets, stride);
            CUDA_CALL(cudaDeviceSynchronize();)
        }


        for(int stride = numBuckets/4; stride > 0; stride /= 2){
            prefix_sum_G1_reverse_second_step<<< num_blocks_prefix_sum_G1, threads_per_block>>>(buckets, resultGPU, numBuckets, stride, true);
            CUDA_CALL(cudaDeviceSynchronize();)
        }


        pippengerMSMG1_unit2_final_add <<<1, 128>>>(resultGPU, buckets);
        CUDA_CALL(cudaDeviceSynchronize());


        if(k > 0){
            pippengerMSMG1_unit3 <<<1, 128>>>(resultGPU, c);
            CUDA_CALL(cudaDeviceSynchronize());
        }
        CUDA_CALL(cudaFree(buckets));
        CUDA_CALL(cudaFree(bucketCounter));
        CUDA_CALL(cudaFree(bucketIndex));
        CUDA_CALL(cudaFree(bucketsBeforeAggregate));

    }

    CUDA_CALL( cudaMemcpy(outputArray, resultGPU,  sizeof(BN254G1), cudaMemcpyDeviceToHost); )
    CUDA_CALL(cudaDeviceSynchronize());
    cudaError_t error = cudaGetLastError();
    if(error != cudaSuccess)
    {
        printf("CUDA error: %s\n", cudaGetErrorString(error));
        exit(-1);
    }
    //this outputArray should only contain one BN254G1 element.

    CUDA_CALL(cudaFree(inputScalarArrayGPU));
    CUDA_CALL(cudaFree(inputBaseArrayGPU));

}




/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseSerialMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseSerialMSMNativeHelper
  (JNIEnv * env, jclass obj,  jbyteArray multiplesOfBaseXYZ, jbyteArray bigScalarsArrayInput, jint batch_size){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");


    vector<Scalar> bigScalarArray = vector<Scalar>(batch_size, Scalar());
    vector<BN254G1> multiplesOfBasePtrArray = vector<BN254G1>(batch_size, BN254G1());


    int len = 32; //254-bit BN254.

    char* baseByteArrayXYZ = (char*)env->GetByteArrayElements(multiplesOfBaseXYZ, NULL);
    for(int i =0; i < batch_size; i++){
        char* tmp = (char*)multiplesOfBasePtrArray[i].X._limbs;
        memcpy(tmp, &baseByteArrayXYZ[i * 3 * len],len);
        char* tmp2 = (char*)multiplesOfBasePtrArray[i].Y._limbs;
        memcpy(tmp2,  &baseByteArrayXYZ[(i *3 +1) * len],len);
        char* tmp3 = (char*)multiplesOfBasePtrArray[i].Z._limbs;
        memcpy(tmp3, &baseByteArrayXYZ[(i*3 +2) * len],len);
    }


    char* scalarBytes = (char*)env->GetByteArrayElements(bigScalarsArrayInput, NULL);
    for(int i =0; i < batch_size; i++){
        char* tmp = (char*)&bigScalarArray[i]._limbs;
        memcpy(tmp, &scalarBytes[i * len ], len);
    }


    //lianke: we know in advance that the numBits will be at most 254. we hard encode it.

    BN254G1* resultCPU = new BN254G1[1];
    memset(resultCPU, 0, sizeof(BN254G1));
    pippengerMSMG1(bigScalarArray, multiplesOfBasePtrArray, resultCPU);

    jbyteArray resultByteArray = env->NewByteArray((jsize)sizeof(BN254G1));

    env->SetByteArrayRegion(resultByteArray, 0 , sizeof(BN254G1) ,   reinterpret_cast<const jbyte*>(&resultCPU[0]));

    return resultByteArray;

  }






















//old unit testcases
__global__ void test_atomicAdd(int* input, int size){
    const int idx = (blockIdx.x * blockDim.x + threadIdx.x);

    int one = 1;
    atomicAdd((int*)&input[idx % size], one);
}

void test_prefix_sum(int size){
    vector<int> inputCPU(size, 1);
    int* inputGPU;
    CUDA_CALL( cudaMalloc((void**)&inputGPU, sizeof(int) * size);)
    CUDA_CALL(cudaMemcpy(inputGPU, &inputCPU[0], sizeof(int) * size, cudaMemcpyHostToDevice);)
    int threads_per_block = 128;
    int num_blocks_prefix_sum = (size + threads_per_block - 1)/threads_per_block;
    //test_atomicAdd<<<num_blocks_prefix_sum * 16, threads_per_block>>> (inputGPU, size);
    for(int stride = 1; stride <= size/2; stride = stride *2){
        prefix_sum_first_step<<< num_blocks_prefix_sum, threads_per_block>>>(inputGPU, size, stride);
        CUDA_CALL(cudaDeviceSynchronize();)
    }
    for(int stride = size/4; stride > 0; stride /= 2){
        prefix_sum_second_step<<< num_blocks_prefix_sum, threads_per_block>>>(inputGPU, size, stride);
        CUDA_CALL(cudaDeviceSynchronize();)
    }
    vector_print<<<1, 128>>> (inputGPU, size);
    CUDA_CALL(cudaDeviceSynchronize();)

}