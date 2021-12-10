#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <vector>
#include <string>
#include <bitset>
#define CHECK_BIT(var,pos) (((var)>>(pos)) & 1)

using namespace std;

class BigInt{
public:
    char* bytes; //we should have 32 bytes
    static const int capacity = 32;
    int len = 0; // number of bytes
    //TODO add modulus here
    BigInt();

    BigInt(char *, int len);
    BigInt(int val);
    //BigInt(const BigInt& val);

    //Direct assignment
    BigInt &operator=(const BigInt &);
 
    //Post/Pre - Incrementation
    BigInt &operator++();
    BigInt operator++(int temp);
    BigInt &operator--();
    BigInt operator--(int temp);
 
    //Addition and Subtraction
    friend BigInt &operator+=(BigInt &, const BigInt &);
    friend BigInt operator+(const BigInt &, const BigInt &);
    friend BigInt operator-(const BigInt &, const BigInt &);
    friend BigInt &operator-=(BigInt &, const BigInt &);
 
    //Comparison operators
    friend bool operator==(const BigInt &, const BigInt &);
    friend bool operator!=(const BigInt &, const BigInt &);
 
    friend bool operator>(const BigInt &, const BigInt &);
    friend bool operator>=(const BigInt &, const BigInt &);
    friend bool operator<(const BigInt &, const BigInt &);
    friend bool operator<=(const BigInt &, const BigInt &);

    //Multiplication and Division
    friend BigInt &operator*=(BigInt &, const BigInt &);
    friend BigInt operator*(const BigInt &, const BigInt &);

    //Modulo
    friend BigInt operator%(const BigInt &, const BigInt &);
    friend BigInt &operator%=(BigInt &, const BigInt &);

    //Power Function
    friend BigInt &operator^=(BigInt &,const BigInt &);
    friend BigInt operator^(BigInt &, const BigInt &);

    //Read and Write
    void print();
    void printAddress();

    bool isZero();
    bool isOne();
    int bitLength();
    bool testBit(int index); //same with the java BigInteger testBit
};

BigInt::BigInt(){
    bytes = new char[capacity];
    memset(bytes, 0, capacity);
}


BigInt::BigInt(int val){
    bytes = new char[capacity];
    memcpy(bytes, &val, capacity - sizeof(int));
}

// BigInt::BigInt(const BigInt& val){
//     bytes = new char[capacity];
//     memcpy(bytes, &val.bytes, capacity);
// }

BigInt::BigInt(char *s, int len){
    bytes = s;
    len = len; 
}

int BigInt::bitLength(){

}

bool BigInt::isZero(){
    //TODO lianke test its correctness
    char testblock[BigInt::capacity];
    memset(testblock, 0, sizeof(testblock));
    return memcmp(testblock, bytes, BigInt::capacity) == 0;
}

bool BigInt::isOne(){
    //TODO lianke test its correctness
    BigInt one(1);
    return *this == one;
}


void BigInt::print(){
    for (int i = 0; i < capacity; i++){
        printf("%hhx |", bytes[i]);
    }
    printf("\n");
    return ;
}

void BigInt::printAddress(){
    printf("0%x\n", bytes);
    return ;
}

bool BigInt::testBit(int index){
    int byte_index = 31 - index / 8; //TODO check its correctness
    int byte_offset = index % 8;
    //printf("%hhx %hhx ", bytes[byte_index], byte_offset);
    return CHECK_BIT(bytes[byte_index], byte_offset);
}

BigInt &BigInt::operator=(const BigInt &a){
    bytes = a.bytes;
    return *this;
}

bool operator==(const BigInt &a, const BigInt &b){
    return memcmp(a.bytes, b.bytes, BigInt::capacity) == 0;
}

BigInt operator+(BigInt &a,const BigInt& b){

    BigInt tmp;
    char mask = 1 << 7;
    bool carry1 = false;
    bool carry2 = false;
    char one = 1;
    char zero = 0;
    for(int i = BigInt::capacity - 1; i > 0; i--){
        tmp.bytes[i] = a.bytes[i] + b.bytes[i] + (carry1||carry2);
        carry1 = ((a.bytes[i] & mask) == mask) && ((b.bytes[i] & mask) == mask);
        carry2 = (((a.bytes[i] & mask) == mask) || ((b.bytes[i] & mask) == mask)) && ((tmp.bytes[i] & mask) !=mask);

        // printf("\n\n  a : %hhx", a.bytes[i]);
        // printf("   b : %hhx", b.bytes[i]);
        // printf("   output : %hhx carry1 %d, carry2 %d \n\n", tmp.bytes[i], carry1, carry2);

    }
    return tmp;
}

//TODO lianke implement several other functions