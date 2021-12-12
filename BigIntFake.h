
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



class BigIntFake{
public:
    uint32_t bytes[16]; //we should have 64 bytes
    static const int capacity = 16;
    static const int num_of_bytes = 64;
    static const int bits_per_word = 32;

    static const long LONG_MASK = 0xffffffffL;

    int len = 0; // number of bytes
    //TODO add modulus here
    BigIntFake();
    BigIntFake(string val);
    BigIntFake(string val, int radix);
    BigIntFake(uint32_t val);
    BigIntFake(int val);

    BigIntFake(unsigned long long val);
    //BigIntFake(const BigIntFake& val);
    
    //Direct assignment
    BigIntFake &operator=(const BigIntFake &);
 
    //Post/Pre - Incrementation
    BigIntFake &operator++();
    BigIntFake operator++(int temp);
    BigIntFake &operator--();
    BigIntFake operator--(int temp);
 
    //Addition and Subtraction
    friend BigIntFake &operator+=(BigIntFake &, const BigIntFake &);
    friend BigIntFake operator+(const BigIntFake &, const BigIntFake &);
    friend BigIntFake operator-(const BigIntFake &, const BigIntFake &);
    friend BigIntFake &operator-=(BigIntFake &, const BigIntFake &);
 
    //Comparison operators
    friend bool operator==(const BigIntFake &, const BigIntFake &);
    friend bool operator!=(const BigIntFake &, const BigIntFake &);
 
    friend bool operator>(const BigIntFake &, const BigIntFake &);
    friend bool operator>=(const BigIntFake &, const BigIntFake &);
    friend bool operator<(const BigIntFake &, const BigIntFake &);
    friend bool operator<=(const BigIntFake &, const BigIntFake &);

    //Multiplication and Division
    friend BigIntFake &operator*=(BigIntFake &, const BigIntFake &);
    friend BigIntFake operator*(const BigIntFake &, const BigIntFake &);

    //Modulo
    friend BigIntFake operator%(const BigIntFake &, const BigIntFake &);
    friend BigIntFake &operator%=(BigIntFake &, const BigIntFake &);

    //Power Function
    friend BigIntFake &operator^=(BigIntFake &,const BigIntFake &);
    friend BigIntFake operator^(BigIntFake &, const BigIntFake &);

    //Read and Write
    void printBinary();
    void printDigits();
    void printAddress();

    void printHex();

    bool isZero();
    bool isOne();
    int bitLength();
    bool testBit(int index); //same with the java BigInteger testBit
    BigIntFake mod(BigIntFake modulus);
};

//const BigIntFake FqModulusParameter = BigIntFake("1532495540865888858358347027150309183618765510462668801");

BigIntFake::BigIntFake(){
    memset(bytes, 0, num_of_bytes);
}


//TODO lianke implement modulus 

BigIntFake::BigIntFake(uint32_t val){
    memcpy(bytes, &val, num_of_bytes - sizeof(uint32_t));
}

BigIntFake::BigIntFake(int val)
{
    memcpy(bytes, &val, num_of_bytes - sizeof(int));
}

BigIntFake::BigIntFake(unsigned long long val){
    memcpy(bytes, &val, num_of_bytes - sizeof(unsigned long long ));
}



bool BigIntFake::isZero(){
    //TODO lianke test its correctness
    uint32_t testblock[BigIntFake::capacity];
    memset(testblock, 0, sizeof(testblock));
    return memcmp(testblock, bytes, num_of_bytes) == 0;
}

bool BigIntFake::isOne(){
    //TODO lianke test its correctness
    BigIntFake one(1);
    return *this == one;
}


void BigIntFake::printBinary(){
    for (int i = 0; i < capacity; i++){
        std::bitset<32> tmp(bytes[i]);
        cout << tmp;
    }
    printf("\n");
    return ;
}

void BigIntFake::printHex(){
    for (int i = 0; i < capacity; i++){
        printf("%X", bytes[i]);

    }
    printf("\n");
    return ;
}

BigIntFake::BigIntFake(string val) : BigIntFake(val, 10){
}

void BigIntFake::printAddress(){
    printf("0%x\n", bytes);
    return ;
}

bool BigIntFake::testBit(int n){
    //TODO lianke this is wrong now, due to we change the internal storage type.
    int byte_index = BigIntFake::capacity - 1 - (n / BigIntFake::bits_per_word); 
    int byte_offset =  n % BigIntFake::bits_per_word;
    //printf("%hhx %hhx ", bytes[byte_index], byte_offset);
    
    return CHECK_BIT(bytes[byte_index], byte_offset);
    //cout << ((bytes[n >> 5] & (1 << (n & 31))) != 0) <<endl;
    //return (bytes[n >> 5] & (1 << (n & 31))) != 0;
}

BigIntFake &BigIntFake::operator=(const BigIntFake &a){
    memcpy(bytes, a.bytes, BigIntFake::num_of_bytes);
    return *this;
}

bool operator==(const BigIntFake &a, const BigIntFake &b){
    return memcmp(a.bytes, b.bytes, BigIntFake::num_of_bytes) == 0;
}

// BigIntFake operator+(BigIntFake &a,const BigIntFake& b){

//     BigIntFake tmp;
//     uint32_t mask = 1 << 31;
//     bool carry1 = false;
//     bool carry2 = false;
//     char one = 1;
//     char zero = 0;

//     for(int i = BigIntFake::capacity - 1; i >= 0; i--){
//         tmp.bytes[i] = a.bytes[i] + b.bytes[i] + (carry1||carry2);
//         carry1 = ((a.bytes[i] & mask) == mask) && ((b.bytes[i] & mask) == mask);
//         carry2 = (((a.bytes[i] & mask) == mask) || ((b.bytes[i] & mask) == mask)) && ((tmp.bytes[i] & mask) !=mask);

//         // printf("\n\n  a : %hhx", a.bytes[i]);
//         // printf("   b : %hhx", b.bytes[i]);
//         // printf("   output : %hhx carry1 %d, carry2 %d \n\n", tmp.bytes[i], carry1, carry2);

//     }
//     return tmp;
// }



BigIntFake operator+(BigIntFake &a, const BigIntFake& b) {
    BigIntFake result;
    uint64_t temp = 0;
    bool carry = false;
    //lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for(int i = BigIntFake::capacity - 1; i >= 0; i--) {
        //cout << a.bytes[i] << " " << b.bytes[i] << " " <<temp << endl;
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        result.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigIntFake::bits_per_word != 0);
    }
    return result;
}



// BigIntFake &operator+=(BigIntFake & a, const BigIntFake & b){
//     uint64_t temp = 0;
//     bool carry = false;
//     //lianke: we only use the lower half. the higher half is for storing larger multiplication results.
//     for(int i = BigIntFake::capacity/2 - 1; i >= 0; i--) {
//         temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
//         a.bytes[i] = (uint32_t)temp;
//         carry = (temp >> 32 != 0);
//     }
//     return a;
// }

// BigIntFake BigIntFake::mod(BigIntFake modulus){
    
// }


// BigIntFake operator*(BigIntFake &a,const BigIntFake& b){

//     BigIntFake result;
//     //TODO lianke implement Karatsuba multiplication.
//     uint16_t temp[BigIntFake::capacity] = {0}; 
//     for(int i = BigIntFake::capacity - 1; i >= BigIntFake::capacity/2; i--){
//         for(int j = BigIntFake::capacity - 1; j>= BigIntFake::capacity/2; j--){
//             temp[i + j - BigIntFake::capacity + 1] += (uint16_t)a.bytes[i] * b.bytes[j];
//         }
//     }
//     uint16_t tmp = 0;
//     uint16_t carry = 0;
//     for(int i = BigIntFake::capacity - 1; i >= 0; i--){
//         tmp = temp[i] + carry;
//         result.bytes[i] = (char)tmp;
//         carry = (tmp >> 32);
//     }
//     return result;
// }





// // Multiply x array times word y in place, and add word z
// void destructiveMulAdd(uint32_t x[], int y, int z) {
//     // Perform the multiplication word by word
//     long ylong = y & BigIntFake::LONG_MASK;
//     long zlong = z & BigIntFake::LONG_MASK;
//     int len = BigIntFake::capacity;

//     unsigned long long product = 0;
//     unsigned long carry = 0;
//     for (int i = len-1; i >= 0; i--) {
//         product = ylong * (x[i] & BigIntFake::LONG_MASK) + carry;
//         x[i] = (uint32_t)product;
//         carry = product >> 32;
//     }

//     // Perform the addition
//     unsigned long sum = (x[len-1] & BigIntFake::LONG_MASK) + zlong;
//     x[len-1] = (uint32_t)sum;
//     carry = sum >> 32;
//     for (int i = len-2; i >= 0; i--) {
//         sum = (x[i] & BigIntFake::LONG_MASK) + carry;
//         x[i] = (uint32_t)sum;
//         carry = sum >> 32;
//     }
// }



// BigIntFake::BigIntFake(string val, int radix) {
//     int cursor = 0;
//     len = val.size();
//     memset(bytes, 0, BigIntFake::num_of_bytes);
//     //to simplify implementation, we do not do input checking here.

//     // Process first (potentially short) digit group
//     int firstGroupLen = len % digitsPerInt[radix];
//     if (firstGroupLen == 0)
//         firstGroupLen = digitsPerInt[radix];
//     string group = val.substr(cursor, firstGroupLen);
//     cursor+= firstGroupLen;
//     bytes[BigIntFake::capacity - 1] = stoi(group);
//     //cout << group << endl;

//     // Process remaining digit groups
//     int superRadix = intRadix[radix];
//     int groupVal = 0;
//     while (cursor < len) {
//         //cout << cursor << " " << cursor + digitsPerInt[radix] << endl;
//         group = val.substr(cursor, digitsPerInt[radix]);
//         cursor += digitsPerInt[radix];
//         //cout << "read in :" << group << endl;
//         groupVal = stoi(group);
//         destructiveMulAdd(bytes, superRadix, groupVal);
//     }
// }

