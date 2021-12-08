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

    bool testBit(int index); //same with the java BigInteger testBit
};

BigInt::BigInt(){
    bytes = new char[capacity];
    memset(bytes, 0, capacity);
}



BigInt::BigInt(char *s, int len){
    bytes = s;
    len = len; 
}


void BigInt::print(){
    for (int i = 0; i < capacity; i++){
        printf("%hhx |", bytes[i]);
    }
    printf("\n");
    return ;
}

bool BigInt::testBit(int index){
    int byte_index = index / 8; //TODO might be reverse
    int byte_offset = index % 8;
    return CHECK_BIT(bytes[byte_index], byte_offset);
}

BigInt &BigInt::operator=(const BigInt &a){
    bytes = a.bytes;
    return *this;
}

BigInt operator+(BigInt &a,const BigInt& b){

    BigInt tmp;
    bool carry = false;
    for(int i = BigInt::capacity - 1; i > 0; i--){
        tmp.bytes[i] = a.bytes[i] + b.bytes[i];
        //TODO update carry bit
    }
    return tmp;
}

//TODO lianke implement several other functions