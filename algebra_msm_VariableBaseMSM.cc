#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <tuple>
#include <cmath>
#include "algebra_msm_VariableBaseMSM.h"
#include "BigInt.h"




/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseSerialMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
//TODO lianke: I plan to pass the modulus through the JNI func parameter to reduce the overhead.
JNIEXPORT jobject JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseSerialMSMNativeHelper
  (JNIEnv * env, jclass obj, jobject bases, jobject scalars){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");


    jint base_size = env->CallIntMethod(bases, java_util_ArrayList_size);
    jint scalars_size = env->CallIntMethod(scalars, java_util_ArrayList_size);

    vector<BigInt> bigScalarArray = vector<BigInt>(scalars_size, BigInt());
    for(int i =0; i < scalars_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(scalars, java_util_ArrayList_get, i);
        bigScalarArray[i].len = env->GetArrayLength(element);
        memcpy(bigScalarArray[i].bytes + BigInt::num_of_bytes - bigScalarArray[i].len, (char*)env->GetByteArrayElements(element, NULL), bigScalarArray[i].len);
        //bigScalarArray[i].print();
    }

    vector<BigInt> baseArray = vector<BigInt>(base_size, BigInt());
    for(int i =0; i < base_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(bases, java_util_ArrayList_get, i);
        baseArray[i].len = env->GetArrayLength(element);
        memcpy(baseArray[i].bytes + BigInt::num_of_bytes - baseArray[i].len, (char*)env->GetByteArrayElements(element, NULL), baseArray[i].len);
        //baseArray[i].print();
    }

    BigInt acc;//should be init to zero.
    vector<tuple<BigInt, BigInt>> filteredInput;
    int numBits = 0;
    for(int i = 0; i < base_size; i++){
      BigInt scalar = bigScalarArray[i];
      if(scalar.isZero()){
        continue;
      }

      BigInt base = baseArray[i];//TODO lianke this base does not contain modulus yet.
      if(scalar.isOne()){
        acc = acc + base;
      }else{
        filteredInput.push_back(make_tuple(scalar, base));
        numBits = max(numBits, scalar.len);
      }
    }


    if(!filteredInput.empty()){
      //TODO lianke call pippengerMSM here. maybe we can make it as a helper function for further usage.

      int length = filteredInput.size();
      int log2Length =  max(1, (int)log2(length));
      int c = log2Length - (log2Length / 3);
      int numBuckets = 1 << c;
      int numGroups = (numBits + c - 1)/c;
      BigInt zero; //TODO lianke modulus should be std::get<1>(filteredInput[0]) they are fakeG1 or fakeG2.
      vector<BigInt> bucketsModel = vector<BigInt>(numBuckets, zero);
      BigInt result = zero;

      for(int k = numGroups - 1; k >=0; k--){
        if (k < numGroups - 1) {
            for (int i = 0; i < c; i++) {
                result = result + result;
            }

            vector<BigInt> buckets = vector<BigInt>(bucketsModel);

            for (int i = 0; i < length; i++) {
                  int id = 0;
                  for (int j = 0; j < c; j++) {
                      if (std::get<1>(filteredInput[i]).testBit(k * c + j)) {
                          id |= 1 << j;
                      }
                  }

                  if (id == 0) {
                      continue;
                  }

                  // Potentially use mixed addition here.
                  buckets[id] = buckets[id] + std::get<1>(filteredInput[i]);
            }

            BigInt runningSum = zero;
            for(int i = numBuckets - 1; i > 0; i--){
              runningSum = runningSum + buckets[i];
              result = result + runningSum;
            }

        }
      }

      acc = acc + result;
    }
    

  //TODO return acc;

    return scalars;

  }

/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseDoubleMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jobject JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseDoubleMSMNativeHelper
  (JNIEnv * env, jclass obj, jobject base1, jobject base2, jobject scalars){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");


    jint base_size1 = env->CallIntMethod(base1, java_util_ArrayList_size);
    jint base_size2 = env->CallIntMethod(base2, java_util_ArrayList_size);

    jint scalars_size = env->CallIntMethod(scalars, java_util_ArrayList_size);

    // vector<BigInt> bigScalarArray = vector<BigInt>(scalars_size, BigInt());
    // for(int i =0; i < scalars_size; i++){
    //     jbyteArray element = (jbyteArray)env->CallObjectMethod(scalars, java_util_ArrayList_get, i);
    //     bigScalarArray[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    // }

    // vector<BigInt> baseArray1 = vector<BigInt>(base_size1, BigInt());
    // for(int i =0; i < base_size1; i++){
    //     jbyteArray element = (jbyteArray)env->CallObjectMethod(base1, java_util_ArrayList_get, i);
    //     baseArray1[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    // }

    // vector<BigInt> baseArray2 = vector<BigInt>(base_size2, BigInt());
    // for(int i =0; i < base_size2; i++){
    //     jbyteArray element = (jbyteArray)env->CallObjectMethod(base2, java_util_ArrayList_get, i);
    //     baseArray2[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    // }

    // BigInt acc1;//should be init to zero.
    // BigInt acc2;//should be init to zero.

    // vector<tuple<BigInt, BigInt>> filteredInput1;
    // vector<tuple<BigInt, BigInt>> filteredInput2;

    // int numBits = 0;
    // for(int i = 0; i < base_size1; i++){
    //   BigInt scalar = bigScalarArray[i];
    //   if(scalar.isZero()){
    //     continue;
    //   }

    //   BigInt base = baseArray1[i];//TODO lianke this base does not contain modulus yet.
    //   if(scalar.isOne()){
    //     acc1 = acc1 + base;
    //   }else{
    //     filteredInput1.push_back(make_tuple(scalar, base));
    //     numBits = max(numBits, scalar.len);
    //   }
    // }

    // for(int i = 0; i < base_size2; i++){
    //   BigInt scalar = bigScalarArray[i];
    //   if(scalar.isZero()){
    //     continue;
    //   }

    //   BigInt base = baseArray2[i];//TODO lianke this base does not contain modulus yet.
    //   if(scalar.isOne()){
    //     acc2 = acc2 + base;
    //   }else{
    //     filteredInput2.push_back(make_tuple(scalar, base));
    //     numBits = max(numBits, scalar.len);
    //   }
    // }

    // if(!filteredInput1.empty()){
    //   int length = filteredInput1.size();
    //   int log2Length =  max(1, (int)log2(length));
    //   int c = log2Length - (log2Length / 3);
    //   int numBuckets = 1 << c;
    //   int numGroups = (numBits + c - 1)/c;
    //   BigInt zero; //TODO lianke modulus should be std::get<1>(filteredInput[0]) they are fakeG1 or fakeG2.
    //   vector<BigInt> bucketsModel = vector<BigInt>(numBuckets, zero);
    //   BigInt result = zero;
    //   for(int k = numGroups - 1; k >=0; k--){
    //     if (k < numGroups - 1) {
    //         for (int i = 0; i < c; i++) {
    //             result = result + result;
    //         }

    //         vector<BigInt> buckets = vector<BigInt>(bucketsModel);

    //         for (int i = 0; i < length; i++) {
    //               int id = 0;
    //               for (int j = 0; j < c; j++) {
    //                   if (std::get<1>(filteredInput1[i]).testBit(k * c + j)) {
    //                       id |= 1 << j;
    //                   }
    //               }

    //               if (id == 0) {
    //                   continue;
    //               }

    //               // Potentially use mixed addition here.
    //               buckets[id] = buckets[id] + std::get<1>(filteredInput1[i]);
    //         }

    //         BigInt runningSum = zero;
    //         for(int i = numBuckets - 1; i > 0; i--){
    //           runningSum = runningSum + buckets[i];
    //           result = result + runningSum;
    //         }

    //     }
    //   }

    //   acc1 = acc1 + result;
    // }
    
  
    // if(!filteredInput2.empty()){
    //   int length = filteredInput2.size();
    //   int log2Length =  max(1, (int)log2(length));
    //   int c = log2Length - (log2Length / 3);
    //   int numBuckets = 1 << c;
    //   int numGroups = (numBits + c - 1)/c;
    //   BigInt zero; //TODO lianke modulus should be std::get<1>(filteredInput[0]) they are fakeG1 or fakeG2.
    //   vector<BigInt> bucketsModel = vector<BigInt>(numBuckets, zero);
    //   BigInt result = zero;
    //   for(int k = numGroups - 1; k >=0; k--){
    //     if (k < numGroups - 1) {
    //         for (int i = 0; i < c; i++) {
    //             result = result + result;
    //         }

    //         vector<BigInt> buckets = vector<BigInt>(bucketsModel);

    //         for (int i = 0; i < length; i++) {
    //               int id = 0;
    //               for (int j = 0; j < c; j++) {
    //                   if (std::get<1>(filteredInput2[i]).testBit(k * c + j)) {
    //                       id |= 1 << j;
    //                   }
    //               }

    //               if (id == 0) {
    //                   continue;
    //               }

    //               // Potentially use mixed addition here.
    //               buckets[id] = buckets[id] + std::get<1>(filteredInput2[i]);
    //         }

    //         BigInt runningSum = zero;
    //         for(int i = numBuckets - 1; i > 0; i--){
    //           runningSum = runningSum + buckets[i];
    //           result = result + runningSum;
    //         }

    //     }
    //   }

    //   acc2 = acc2 + result;
    // }

    //TODO return acc1 and acc2

    return scalars;
  }