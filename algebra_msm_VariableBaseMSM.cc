#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include "algebra_msm_VariableBaseMSM.h"
#include "BigInteger.h"

/*
 * Class:     algebra_msm_VariableBaseMSM
 * Method:    variableBaseSerialMSMNativeHelper
 * Signature: (Ljava/util/ArrayList;Ljava/util/ArrayList;)Lalgebra/groups/AbstractGroup;
 */
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
        bigScalarArray[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    }

    vector<BigInt> baseArray = vector<BigInt>(base_size, BigInt());
    for(int i =0; i < base_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(bases, java_util_ArrayList_get, i);
        baseArray[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    }


    //TODO lianke do the variable MSM computation.


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

    vector<BigInt> bigScalarArray = vector<BigInt>(scalars_size, BigInt());
    for(int i =0; i < scalars_size; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(scalars, java_util_ArrayList_get, i);
        bigScalarArray[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    }

    vector<BigInt> baseArray1 = vector<BigInt>(base_size1, BigInt());
    for(int i =0; i < base_size1; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(base1, java_util_ArrayList_get, i);
        baseArray1[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    }

    vector<BigInt> baseArray2 = vector<BigInt>(base_size2, BigInt());
    for(int i =0; i < base_size2; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(base2, java_util_ArrayList_get, i);
        baseArray2[i].bytes = (char*)env->GetByteArrayElements(element, NULL);
    }

    //TODO lianke do the variable MSM computation.

    return scalars;
  }