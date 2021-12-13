
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <vector>
#include <string>
#include <bitset>
#define CHECK_BIT(var,pos) (((var)>>(pos)) & 1)

using namespace std;
//TODO we also need to implement BN254, not just for FakeG1 and FakeG2


int digitsPerInt[] = {0, 0, 30, 19, 15, 13, 11,
        11, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6,
        6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5};

int intRadix[] = {0, 0,
        0x40000000, 0x4546b3db, 0x40000000, 0x48c27395, 0x159fd800,
        0x75db9c97, 0x40000000, 0x17179149, 0x3b9aca00, 0xcc6db61,
        0x19a10000, 0x309f1021, 0x57f6c100, 0xa2f1b6f,  0x10000000,
        0x18754571, 0x247dbc80, 0x3547667b, 0x4c4b4000, 0x6b5a6e1d,
        0x6c20a40,  0x8d2d931,  0xb640000,  0xe8d4a51,  0x1269ae40,
        0x17179149, 0x1cb91000, 0x23744899, 0x2b73a840, 0x34e63b41,
        0x40000000, 0x4cfa3cc1, 0x5c13d840, 0x6d91b519, 0x39aa400
    };


class BigInt{
public:
    uint32_t bytes[16]; //we should have 64 bytes
    static const int capacity = 16;
    static const int num_of_bytes = 64;
    static const int bits_per_word = 32;

    static const long LONG_MASK = 0xffffffffL;

    int len = 0; // number of bytes
    //TODO add modulus here
    BigInt();
    BigInt(string val);
    BigInt(string val, int radix);
    BigInt(uint32_t val);
    BigInt(int val);


    BigInt(unsigned long long val);
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
    void printBinary();
    void printDigits();
    void printAddress();

    void printHex();

    static BigInt ZERO();
    static BigInt ONE();
    bool isZero();
    bool isOne();
    int bitLength();
    bool testBit(int index); //same with the java BigInteger testBit
    BigInt mod(BigInt modulus);
};

//const BigInt FqModulusParameter = BigInt("1532495540865888858358347027150309183618765510462668801");

BigInt BigInt::ZERO(){
    return BigInt(0);
}

BigInt BigInt::ONE(){
    return BigInt(1);
}


BigInt::BigInt(){
    memset(bytes, 0, num_of_bytes);
}


//TODO lianke implement modulus 

BigInt::BigInt(uint32_t val){
    memcpy(bytes, &val, num_of_bytes - sizeof(uint32_t));
}

BigInt::BigInt(int val)
{
    memcpy(bytes, &val, num_of_bytes - sizeof(int));
}

BigInt::BigInt(unsigned long long val){
    memcpy(bytes, &val, num_of_bytes - sizeof(unsigned long long ));
}



bool BigInt::isZero(){
    //TODO lianke test its correctness
    uint32_t testblock[BigInt::capacity];
    memset(testblock, 0, sizeof(testblock));
    return memcmp(testblock, bytes, num_of_bytes) == 0;
}

bool BigInt::isOne(){
    //TODO lianke test its correctness
    BigInt one(1);
    return *this == one;
}


void BigInt::printBinary(){
    for (int i = 0; i < capacity; i++){
        std::bitset<32> tmp(bytes[i]);
        cout << tmp;
    }
    printf("\n");
    return ;
}

void BigInt::printHex(){
    for (int i = 0; i < capacity; i++){
        printf("%X", bytes[i]);

    }
    printf("\n");
    return ;
}

BigInt::BigInt(string val) : BigInt(val, 10){

}

void BigInt::printAddress(){
    printf("0%x\n", bytes);
    return ;
}

bool BigInt::testBit(int n){
    //TODO lianke this is wrong now, due to we change the internal storage type.
    int byte_index = BigInt::capacity - 1 - (n / BigInt::bits_per_word); 
    int byte_offset =  n % BigInt::bits_per_word;
    //printf("%hhx %hhx ", bytes[byte_index], byte_offset);
    
    return CHECK_BIT(bytes[byte_index], byte_offset);
    //cout << ((bytes[n >> 5] & (1 << (n & 31))) != 0) <<endl;
    //return (bytes[n >> 5] & (1 << (n & 31))) != 0;
}

BigInt &BigInt::operator=(const BigInt &a){
    memcpy(bytes, a.bytes, BigInt::num_of_bytes);
    return *this;
}

bool operator==(const BigInt &a, const BigInt &b){
    return memcmp(a.bytes, b.bytes, BigInt::num_of_bytes) == 0;
}





BigInt operator+(BigInt &a, const BigInt& b) {
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;
    //lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for(int i = BigInt::capacity - 1; i >= 0; i--) {
        //cout << a.bytes[i] << " " << b.bytes[i] << " " <<temp << endl;
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        result.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigInt::bits_per_word != 0);
    }
    return result;
}



BigInt &operator+=(BigInt & a, const BigInt & b){
    uint64_t temp = 0;
    bool carry = false;
    //lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for(int i = BigInt::capacity/2 - 1; i >= 0; i--) {
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        a.bytes[i] = (uint32_t)temp;
        carry = (temp >> 32 != 0);
    }
    return a;
}

BigInt BigInt::mod(BigInt modulus){
    //TODO lianke

}


BigInt operator*(BigInt &a,const BigInt& b){

    BigInt result;
    //TODO lianke implement Karatsuba multiplication.
    uint16_t temp[BigInt::capacity] = {0}; 
    for(int i = BigInt::capacity - 1; i >= BigInt::capacity/2; i--){
        for(int j = BigInt::capacity - 1; j>= BigInt::capacity/2; j--){
            temp[i + j - BigInt::capacity + 1] += (uint16_t)a.bytes[i] * b.bytes[j];
        }
    }
    uint16_t tmp = 0;
    uint16_t carry = 0;
    for(int i = BigInt::capacity - 1; i >= 0; i--){
        tmp = temp[i] + carry;
        result.bytes[i] = (char)tmp;
        carry = (tmp >> 32);
    }
    return result;
}





// Multiply x array times word y in place, and add word z
void destructiveMulAdd(uint32_t x[], int y, int z) {
    // Perform the multiplication word by word
    long ylong = y & BigInt::LONG_MASK;
    long zlong = z & BigInt::LONG_MASK;
    int len = BigInt::capacity;

    unsigned long long product = 0;
    unsigned long carry = 0;
    for (int i = len-1; i >= 0; i--) {
        product = ylong * (x[i] & BigInt::LONG_MASK) + carry;
        x[i] = (uint32_t)product;
        carry = product >> 32;
    }

    // Perform the addition
    unsigned long sum = (x[len-1] & BigInt::LONG_MASK) + zlong;
    x[len-1] = (uint32_t)sum;
    carry = sum >> 32;
    for (int i = len-2; i >= 0; i--) {
        sum = (x[i] & BigInt::LONG_MASK) + carry;
        x[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
}



BigInt::BigInt(string val, int radix) {
    int cursor = 0;
    len = val.size();
    memset(bytes, 0, BigInt::num_of_bytes);
    //to simplify implementation, we do not do input checking here.

    // Process first (potentially short) digit group
    int firstGroupLen = len % digitsPerInt[radix];
    if (firstGroupLen == 0)
        firstGroupLen = digitsPerInt[radix];
    string group = val.substr(cursor, firstGroupLen);
    cursor+= firstGroupLen;
    bytes[BigInt::capacity - 1] = stoi(group);
    //cout << group << endl;

    // Process remaining digit groups
    int superRadix = intRadix[radix];
    int groupVal = 0;
    while (cursor < len) {
        //cout << cursor << " " << cursor + digitsPerInt[radix] << endl;
        group = val.substr(cursor, digitsPerInt[radix]);
        cursor += digitsPerInt[radix];
        //cout << "read in :" << group << endl;
        groupVal = stoi(group);
        destructiveMulAdd(bytes, superRadix, groupVal);
    }
}

