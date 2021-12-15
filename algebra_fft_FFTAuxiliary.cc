#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <bitset>
#include <vector>
#include <cmath>
#include <chrono>
#include "algebra_fft_FFTAuxiliary.h"
#include "BigInt.h"


int reverseBits(int n, int range) {
    int ans = 0;
    for(int i = range - 1; i >= 0; i--){
        ans |= (n & 1) <<i;
        n>>=1;
    }
    return ans;
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


    vector<BigInt> inputArray = vector<BigInt>(input_len, BigInt());

    for(int i =0; i < input_len; i++){
        jbyteArray element = (jbyteArray)env->CallObjectMethod(inputs, java_util_ArrayList_get, i);
        char* bytes = (char*)env->GetByteArrayElements(element, NULL);
        inputArray[i].len = env->GetArrayLength(element);
        char* tmp = (char*)&inputArray[i].bytes;

        memcpy(tmp +BigInt::num_of_bytes - inputArray[i].len, 
                                    bytes,
                                    inputArray[i].len);
    }

    BigInt omega("0");
    char* bytes = (char*)env->GetByteArrayElements(omegaArray, NULL);
    omega.len = env->GetArrayLength(omegaArray);
    char* tmp = (char*)&omega.bytes;
    memcpy(tmp +BigInt::num_of_bytes - omega.len, 
                            bytes,
                            omega.len);
    

    int logn = log2(input_len);    
    for (uint32_t k = 0; k < input_len; ++k) {
        uint32_t rk = reverseBits(k, logn);

        if (k < rk) {
            BigInt temp = inputArray[k];
            inputArray[k] = inputArray[rk];
            inputArray[rk] = temp;
        }
    }

    int m = 1; // invariant: m = 2^{s-1}
    omega.printBinary();
    //TODO LIANKE i think this part is the key overhead in serialFFT, we should CUDA it.
    for (int s = 1; s <= logn; ++s) {
        // w_m is 2^s-th root of unity now
        cout << "s=" <<s <<endl;
        BigInt w_m = pow(omega, (int)input_len / (2 * m));
        cout <<"finished pow"<<endl;
        w_m.printBinary();
        for (int k = 0; k < input_len; k += 2 * m) {
            BigInt w("1");
            //cout<< "k="<< k <<endl;
            for (int j = 0; j < m; ++j) {
                BigInt t = w * (inputArray[k + j + m]);
                inputArray[k + j + m]= inputArray[k + j] - t;
                inputArray[k + j] = inputArray[k + j] + t;
                w = w * w_m;
            }
        }
        m *= 2;
    }

    jbyteArray resultByteArray = env->NewByteArray((jsize)BigInt::num_of_bytes * input_len);

    for(int i=0; i < input_len;i++){
        env->SetByteArrayRegion(resultByteArray, i * BigInt::num_of_bytes , BigInt::num_of_bytes,   reinterpret_cast<const jbyte*>(inputArray[i].bytes));
    }


    return resultByteArray;

}
