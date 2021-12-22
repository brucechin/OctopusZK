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

// BigInt FqModulusParameter("1532495540865888858358347027150309183618765510462668801");



/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseSerialMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseSerialMSMNativeHelper
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
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        char* tmp = (char*)&bigScalarArray[i].bytes;

        memcpy(tmp +BigInt::num_of_bytes - bigScalarArray[i].len, 
                                bytes,
                                bigScalarArray[i].len);
    }

    vector<BigInt> baseArray = vector<BigInt>(base_size, BigInt());
    for(int i =0; i < base_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(bases, java_util_ArrayList_get, i);
        baseArray[i].len = env->GetArrayLength(element);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        char* tmp = (char*)baseArray[i].bytes;
        memcpy(tmp + BigInt::num_of_bytes - baseArray[i].len, bytes,  baseArray[i].len);
    }


    BigInt acc("0");//should be init to zero.
    // cout << " ----------------acc init-------------"<<endl;
    // acc.printBinary();
    // cout << " ------------------------------------"<<endl;
    vector<tuple<BigInt, BigInt>> filteredInput;
    int numBits = 0;
    for(int i = 0; i < base_size; i++){
      BigInt scalar = bigScalarArray[i];
      if(scalar.isZero()){
        continue;
      }

      BigInt base = baseArray[i];
      if(scalar.isOne()){
        acc = (acc + base);
        acc  %= FqModulusParameter;
      }else{
        filteredInput.push_back(make_tuple(scalar, base));
        numBits = max(numBits, scalar.bitLength());
        // cout << "cpp side add <scalar, base> index:" << i << " \n"; 
        // scalar.printBinary();
        // base.printBinary();
      }

    }

    //cout << "cpp side filteredInput size: " << filteredInput.size() << " numBits : " << numBits << endl;
    // cout << "cpp side current acc is ";
    // acc.printBinary();

    jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes);

    if(!filteredInput.empty()){

      int length = filteredInput.size();
      int log2Length =  max(1, (int)log2(length));
      int c = log2Length - (log2Length / 3);
      int numBuckets = 1 << c;
      int numGroups = (numBits + c - 1)/c;
      BigInt zero("0"); 
      vector<BigInt> bucketsModel = vector<BigInt>(numBuckets, zero);
      BigInt result("0"); 
      //cout << "cpp side length " << length << " log2 length " << log2Length <<" c "  << c <<" numGroups " << numGroups << " numBuckets " << numBuckets << endl;
      for(int k = numGroups - 1; k >=0; k--){
        if (k < numGroups - 1) {
            for (int i = 0; i < c; i++) {
                result = (result + result) ;
                result %= FqModulusParameter; 
            }
        }
      //   if(k > numGroups - 20) { 
      //  cout << "cpp side k=" << k << " after exp result is :";
      //   result.printBinary();
      //   }


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

              // if(k == 80 && i < 10){
              //   cout << "cpp side i=" << i << " after testBit, id is :" << id << endl;
              // }

              // Potentially use mixed addition here.
              buckets[id] = (buckets[id] + std::get<1>(filteredInput[i]));
              buckets[id]  %= FqModulusParameter;
        }
        // if(k == 80){
        //   for(int i = 0; i < numBuckets; i++){
        //     cout <<"cpp side buckets index " << i << ":";
        //     buckets[i].printBinary();
        //   }
        // }

        BigInt runningSum("0");
        for(int i = numBuckets - 1; i > 0; i--){
          // if(k == numGroups - 1){
          //   cout << "cpp side k " << k << " bucket index " << i << " runningSum" << endl;
          //   runningSum.printBinary();
          //   result.printBinary();
          // }
          runningSum = (runningSum + buckets[i]);
          runningSum %= FqModulusParameter;
          result = (result + runningSum);
          result %= FqModulusParameter;
        }
        // if(k < 10){ 
        // cout << "cpp side numGroups:" << numGroups <<" k= " << k << " after adding runningSum result is :";
        // result.printBinary();
        // }
        //env->SetByteArrayRegion(resultByteArray, k * BigInt::num_of_bytes, BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(result.bytes));

      }


      acc = (acc + result);
      acc %= FqModulusParameter;
    }
    // cout << "cpp side final acc" << endl;
    // acc.printBinary();


    env->SetByteArrayRegion(resultByteArray, 0 , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(acc.bytes));

    return resultByteArray;

  }

/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseDoubleMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
JNIEXPORT jbyteArray JNICALL Java_algebra_msm_VariableBaseMSM_variableBaseDoubleMSMNativeHelper
  (JNIEnv * env, jclass obj, jobject base1, jobject base2, jobject scalars){
    jclass java_util_ArrayList      = static_cast<jclass>(env->NewGlobalRef(env->FindClass("java/util/ArrayList")));
    jmethodID java_util_ArrayList_size = env->GetMethodID(java_util_ArrayList, "size", "()I");
    jmethodID java_util_ArrayList_get  = env->GetMethodID(java_util_ArrayList, "get", "(I)Ljava/lang/Object;");


    jint base_size1 = env->CallIntMethod(base1, java_util_ArrayList_size);
    jint base_size2 = env->CallIntMethod(base2, java_util_ArrayList_size);

    jint scalars_size = env->CallIntMethod(scalars, java_util_ArrayList_size);


    vector<BigInt> bigScalarArray = vector<BigInt>(scalars_size, BigInt());
    for(int i =0; i < scalars_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(scalars, java_util_ArrayList_get, i);
        bigScalarArray[i].len = env->GetArrayLength(element);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        char* tmp = (char*)&bigScalarArray[i].bytes;

        memcpy(tmp +BigInt::num_of_bytes - bigScalarArray[i].len, 
                                bytes,
                                bigScalarArray[i].len);
    }

    vector<BigInt> baseArray1 = vector<BigInt>(base_size1, BigInt());
    for(int i =0; i < base_size1; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(base1, java_util_ArrayList_get, i);
        baseArray1[i].len = env->GetArrayLength(element);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        char* tmp = (char*)baseArray1[i].bytes;
        memcpy(tmp + BigInt::num_of_bytes - baseArray1[i].len, bytes,  baseArray1[i].len);
    }


    vector<BigInt> baseArray2 = vector<BigInt>(base_size2, BigInt());
    for(int i =0; i < base_size2; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(base2, java_util_ArrayList_get, i);
        baseArray2[i].len = env->GetArrayLength(element);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        char* tmp = (char*)baseArray2[i].bytes;
        memcpy(tmp + BigInt::num_of_bytes - baseArray2[i].len, bytes,  baseArray2[i].len);
    }



    BigInt acc1("0");
    BigInt acc2("0");

    vector<tuple<BigInt, BigInt>> filteredInput1;
    vector<tuple<BigInt, BigInt>> filteredInput2;

    int numBits = 0;

    for(int i = 0; i < base_size1; i++){
      BigInt scalar = bigScalarArray[i];
      if(scalar.isZero()){
        continue;
      }

      BigInt base = baseArray1[i];
      if(scalar.isOne()){
        acc1 = acc1 + base;
        acc1 %= FqModulusParameter;

      }else{
        filteredInput1.push_back(make_tuple(scalar, base));
        numBits = max(numBits, scalar.bitLength());
      }
    }

    for(int i = 0; i < base_size2; i++){
      BigInt scalar = bigScalarArray[i];
      if(scalar.isZero()){
        continue;
      }

      BigInt base = baseArray2[i];
      if(scalar.isOne()){
        acc2 = acc2 + base;
        acc2 %= FqModulusParameter;

      }else{
        filteredInput2.push_back(make_tuple(scalar, base));
        //numBits = max(numBits, scalar.bitLength());
      }
    }

    // for(int i = 0; i < base_size1; i++){
    //   BigInt scalar = bigScalarArray[i];
    //   if(scalar.isZero()){
    //     continue;
    //   }

    //   BigInt base1 = baseArray1[i];
    //   BigInt base2 = baseArray2[i];
    //   if(scalar.isOne()){
    //     acc1 = acc1 + base1;
    //     acc1 %= FqModulusParameter;
    //     acc2 = acc2 + base2;
    //     acc2 %= FqModulusParameter;
    //   }else{
    //     filteredInput1.push_back(make_tuple(scalar, base1));
    //     filteredInput2.push_back(make_tuple(scalar, base2));
    //     numBits = max(numBits, scalar.bitLength());
    //   }
    // }


    if(!filteredInput1.empty()){
      int length = filteredInput1.size();
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
                result = result + result;
                result %= FqModulusParameter; 
            }
        }
        vector<BigInt> buckets = vector<BigInt>(bucketsModel);

        for (int i = 0; i < length; i++) {
              int id = 0;
              for (int j = 0; j < c; j++) {
                  if (std::get<0>(filteredInput1[i]).testBit(k * c + j)) {
                      id |= 1 << j;
                  }
              }

              if (id == 0) {
                  continue;
              }

              // Potentially use mixed addition here.
              buckets[id] = buckets[id] + std::get<1>(filteredInput1[i]);
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

      acc1 = acc1 + result;
      acc1 %= FqModulusParameter;

    }
    
  
    if(!filteredInput2.empty()){
      int length = filteredInput2.size();
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
                result = result + result;
                result %= FqModulusParameter; 
            }
        }

            vector<BigInt> buckets = vector<BigInt>(bucketsModel);

            for (int i = 0; i < length; i++) {
                  int id = 0;
                  for (int j = 0; j < c; j++) {
                      if (std::get<0>(filteredInput2[i]).testBit(k * c + j)) {
                          id |= 1 << j;
                      }
                  }

                  if (id == 0) {
                      continue;
                  }

                  // Potentially use mixed addition here.
                  buckets[id] = buckets[id] + std::get<1>(filteredInput2[i]);
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

      acc2 = acc2 + result;
      acc2 %= FqModulusParameter;

    }

    // cout << "cpp side final acc1 and acc2"<<endl;
    // acc1.printBinary();
    // acc2.printBinary();

    jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes * 2 );
    env->SetByteArrayRegion(resultByteArray, 0 , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(acc1.bytes));
    env->SetByteArrayRegion(resultByteArray, BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(acc2.bytes));



    return resultByteArray;
  }