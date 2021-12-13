#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <chrono>
#include "algebra_msm_FixedBaseMSM.h"
//#include "BigInteger.h"
#include "BigInt.h"

//TODO lianke: G1 and G2 MSM window table generation can be moved to cpp side too.



using namespace std;
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_FixedBaseMSM_batchMSMNativeHelper
  (JNIEnv *env, jclass obj, jint outerc, jint windowSize, jobject multiplesOfBase, jobject bigScalars)
{

  jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
  jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
  jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");

  jint out_len = env->CallIntMethod(multiplesOfBase, java_util_ArrayList_size);
  jint inner_len = env->CallIntMethod(env->CallObjectMethod(multiplesOfBase, java_util_ArrayList_get, 0), java_util_ArrayList_size);
  
  jint batch_size = env->CallIntMethod(bigScalars, java_util_ArrayList_size);
  //cout << "cpp side batch size: " << batch_size << endl;


  auto start = std::chrono::steady_clock::now();
  vector<BigInt> bigScalarArray = vector<BigInt>(batch_size, BigInt());

  vector<vector<BigInt>> multiplesOfBasePtrArray = vector<vector<BigInt>>(out_len, vector<BigInt>(inner_len, BigInt()));

    auto end = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "BigInt allocation elapsed time: " << elapsed_seconds.count() << "s\n";



  start = std::chrono::steady_clock::now();
  for(int i =0; i < batch_size; i++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      bigScalarArray[i].len = env->GetArrayLength(element);
      //bigScalarArrayFake[i].len = env->GetArrayLength(element);
      char* tmp = (char*)&bigScalarArray[i].bytes;

      memcpy(tmp +BigInt::num_of_bytes - bigScalarArray[i].len, 
                                bytes,
                                bigScalarArray[i].len);


      // cout << i << " " <<bigScalarArray[i].len<<endl;
      // for (int j = 0; j < bigScalarArray[i].len ; j++){
      //   std::bitset<8> tmp(bytes[j]);
      //   cout << tmp;
      // }
      // cout <<endl;
      // bigScalarArray[i].printBinary();
      // bigScalarArrayFake[i].printBinary();
      // cout << "-----------------------------"<<endl;
      //}


  }

    //TODO lianke parallelize it
  //TODO lianke use GetIntArrayElements instead of GetByteArrayElements
  for(int i = 0; i < out_len;i++){
    for(int j = 0; j < inner_len; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      multiplesOfBasePtrArray[i][j].len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray[i][j].bytes;
      memcpy(tmp + BigInt::num_of_bytes - multiplesOfBasePtrArray[i][j].len, bytes,  multiplesOfBasePtrArray[i][j].len);


      //multiplesOfBasePtrArrayFake[i][j].printBinary();
      // cout << "--------multiplesOfBase cmp------------"<<endl;
                      
      // cout << memcmp(multiplesOfBasePtrArray[i][j].bytes, multiplesOfBasePtrArrayFake[i][j].bytes, BigInt::num_of_bytes) <<endl;

      //cout << i << " " <<j  << " bytes len:" <<multiplesOfBasePtrArray[i][j].len <<endl;
      //multiplesOfBasePtrArray[i][j].print();
    }
  }
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end-start;
    std::cout << "Read from JVM elapsed time: " << elapsed_seconds.count() << "s\n";


  start = std::chrono::steady_clock::now();
  jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes * batch_size);
  for(int batch_index = 0; batch_index < batch_size; batch_index++){
    BigInt res = multiplesOfBasePtrArray[0][0];//TODO lianke this assignment has a problem. 

    for (int outer = 0; outer < outerc; ++outer) {
        int inner = 0;

        for (int i = 0; i < windowSize; ++i) {
            if (bigScalarArray[batch_index].testBit(outer * windowSize + i)) { //Returns true if and only if the designated bit is set.
                inner |= 1 << i;
            }
          //   if(bigScalarArrayFake[batch_index].testBit(outer * windowSize + i) != bigScalarArray[batch_index].testBit(outer * windowSize + i)){
          //     cout << "testBit error " <<batch_index << " " <<  outer * windowSize + i << endl;
          //     cout << bigScalarArrayFake[batch_index].testBit(outer * windowSize + i) << " " << bigScalarArray[batch_index].testBit(outer * windowSize + i) <<endl;
          // }
        }


        res = res + multiplesOfBasePtrArray[outer][inner];

    }  





    env->SetByteArrayRegion(resultByteArray, batch_index * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(res.bytes));
  }
    end = std::chrono::steady_clock::now();
    elapsed_seconds = end-start;
    std::cout << "C++ Compute elapsed time: " << elapsed_seconds.count() << "s\n";
  

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


  auto start = std::chrono::steady_clock::now();
  vector<BigInt> bigScalarArray = vector<BigInt>(batch_size, BigInt());
  vector<vector<BigInt>> multiplesOfBasePtrArray1 = vector<vector<BigInt>>(out_len1, vector<BigInt>(inner_len1, BigInt()));
  vector<vector<BigInt>> multiplesOfBasePtrArray2 = vector<vector<BigInt>>(out_len2, vector<BigInt>(inner_len2, BigInt()));

  auto end = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = end-start;
  std::cout << "doubleBatchMSM BigInt allocation elapsed time: " << elapsed_seconds.count() << "s\n";




  start = std::chrono::steady_clock::now();
  for(int i =0; i < batch_size; i++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(bigScalars, java_util_ArrayList_get, i);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      bigScalarArray[i].len = env->GetArrayLength(element);
      char* tmp = (char*)&bigScalarArray[i].bytes;

      memcpy(tmp +BigInt::num_of_bytes - bigScalarArray[i].len, 
                                bytes,
                                bigScalarArray[i].len);

  }



  //TODO parallelize it
  for(int i = 0; i < out_len1;i++){
    for(int j = 0; j < inner_len1; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase1, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      multiplesOfBasePtrArray1[i][j].len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray1[i][j].bytes;
      memcpy(tmp + BigInt::num_of_bytes - multiplesOfBasePtrArray1[i][j].len, bytes,  multiplesOfBasePtrArray1[i][j].len);
    }
  }

  for(int i = 0; i < out_len2;i++){
    for(int j = 0; j < inner_len2; j++){
      jbyteArray element = (jbyteArray)env->CallObjectMethod(env->CallObjectMethod(multiplesOfBase2, java_util_ArrayList_get, i), java_util_ArrayList_get, j);
      char* bytes = (char*)env->GetByteArrayElements(element, NULL);
      multiplesOfBasePtrArray2[i][j].len = env->GetArrayLength(element);
      char* tmp = (char*)multiplesOfBasePtrArray2[i][j].bytes;
      memcpy(tmp + BigInt::num_of_bytes - multiplesOfBasePtrArray2[i][j].len, bytes,  multiplesOfBasePtrArray2[i][j].len);    }
  }

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end-start;
    std::cout << "doubleBatchMSM Read from JVM elapsed time: " << elapsed_seconds.count() << "s\n";


  start = std::chrono::steady_clock::now();
  jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes * batch_size * 2);
  for(int batch_index = 0; batch_index < batch_size; batch_index++){
    BigInt res1 = multiplesOfBasePtrArray1[0][0];
    for (int outer = 0; outer < outerc1; ++outer) {
        int inner = 0;
        for (int i = 0; i < windowSize1; ++i) {
            if (bigScalarArray[batch_index].testBit(outer * windowSize1 + i)) { //Returns true if and only if the designated bit is set.
                inner |= 1 << i;
            }
        }
        res1 = res1 + multiplesOfBasePtrArray1[outer][inner];
    }



    BigInt res2 = multiplesOfBasePtrArray2[0][0];
    for (int outer = 0; outer < outerc2; ++outer) {
        int inner = 0;
        for (int i = 0; i < windowSize2; ++i) {
            if (bigScalarArray[batch_index].testBit(outer * windowSize2 + i)) { //Returns true if and only if the designated bit is set.
                inner |= 1 << i;
            }
        }
        res2 = res2 + multiplesOfBasePtrArray2[outer][inner];
    }

    env->SetByteArrayRegion(resultByteArray, 2 * batch_index * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(res1.bytes));
    env->SetByteArrayRegion(resultByteArray, (2 * batch_index + 1) * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(res2.bytes));

  }

    end = std::chrono::steady_clock::now();
    elapsed_seconds = end-start;
    std::cout << "doubleBatchMSM C++ Compute elapsed time: " << elapsed_seconds.count() << "s\n";
  

  return resultByteArray;

  }




