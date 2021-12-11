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
    uint32_t* bytes; //we should have 32 bytes
    static const int capacity = 8;
    static const int num_of_bytes = 32;

    int len = 0; // number of bytes
    //TODO add modulus here
    BigInt();
    BigInt(string val);
    BigInt(unsigned int val);
    BigInt(unsigned long long val);
    BigInt(BigInt val);
    
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
const BigInt FqModulusParameter = BigInt("1532495540865888858358347027150309183618765510462668801");

BigInt::BigInt(){
    bytes = new uint32_t[BigInt::capacity];
    memset(bytes, 0, num_of_bytes);
}


//TODO lianke implement modulus 

BigInt::BigInt(unsigned int val){
    bytes = new uint32_t[BigInt::capacity];
    memcpy(bytes, &val, num_of_bytes - sizeof(unsigned int));
}

BigInt::BigInt(unsigned long long val){
    bytes = new uint32_t[BigInt::capacity];
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


void BigInt::print(){
    for (int i = 0; i < capacity; i++){
        std::bitset<32> tmp(bytes[i]);
        cout << tmp << "|";
    }
    printf("\n");
    return ;
}

void BigInt::printAddress(){
    printf("0%x\n", bytes);
    return ;
}

bool BigInt::testBit(int n){
    //TODO lianke this is wrong now, due to we change the internal storage type.
    // int byte_index = BigInt::capacity - 1 - index / 8; 
    // int byte_offset = index % 8;
    // //printf("%hhx %hhx ", bytes[byte_index], byte_offset);
    // return CHECK_BIT(bytes[byte_index], byte_offset);
    return (bytes[n >>> 5] & (1 << (n & 31))) != 0;
}

BigInt &BigInt::operator=(const BigInt &a){
    memcpy(bytes, a.bytes, BigInt::num_of_bytes);
    return *this;
}

bool operator==(const BigInt &a, const BigInt &b){
    return memcmp(a.bytes, b.bytes, BigInt::num_of_bytes) == 0;
}

// BigInt operator+(BigInt &a,const BigInt& b){

//     BigInt tmp;
//     char mask = 1 << 7;
//     bool carry1 = false;
//     bool carry2 = false;
//     char one = 1;
//     char zero = 0;
//     int max_len = max(a.len, b.len);
//     for(int i = BigInt::capacity - 1; i > 0; i--){
//         tmp.bytes[i] = a.bytes[i] + b.bytes[i] + (carry1||carry2);
//         carry1 = ((a.bytes[i] & mask) == mask) && ((b.bytes[i] & mask) == mask);
//         carry2 = (((a.bytes[i] & mask) == mask) || ((b.bytes[i] & mask) == mask)) && ((tmp.bytes[i] & mask) !=mask);

//         // printf("\n\n  a : %hhx", a.bytes[i]);
//         // printf("   b : %hhx", b.bytes[i]);
//         // printf("   output : %hhx carry1 %d, carry2 %d \n\n", tmp.bytes[i], carry1, carry2);

//     }
//     return tmp;
// }

BigInt &operator+=(BigInt & lvalue, const BigInt & rvalue){

}



BigInt operator*(BigInt &a,const BigInt& b){

    BigInt tmp;

    return tmp;
}

//TODO lianke implement several other functions

// Multiply x array times word y in place, and add word z
    void destructiveMulAdd(int[] x, int y, int z) {
    // Perform the multiplication word by word
    long ylong = y & LONG_MASK;
    long zlong = z & LONG_MASK;
    int len = x.length;

    long product = 0;
    long carry = 0;
    for (int i = len-1; i >= 0; i--) {
        product = ylong * (x[i] & LONG_MASK) + carry;
        x[i] = (int)product;
        carry = product >>> 32;
    }

    // Perform the addition
    long sum = (x[len-1] & LONG_MASK) + zlong;
    x[len-1] = (int)sum;
    carry = sum >>> 32;
    for (int i = len-2; i >= 0; i--) {
        sum = (x[i] & LONG_MASK) + carry;
        x[i] = (int)sum;
        carry = sum >>> 32;
    }
}


/**
 * Translates the String representation of a BigInteger in the
 * specified radix into a BigInteger.  The String representation
 * consists of an optional minus or plus sign followed by a
 * sequence of one or more digits in the specified radix.  The
 * character-to-digit mapping is provided by {@code
 * Character.digit}.  The String may not contain any extraneous
 * characters (whitespace, for example).
 *
 * @param val String representation of BigInteger.
 * @param radix radix to be used in interpreting {@code val}.
 * @throws NumberFormatException {@code val} is not a valid representation
 *         of a BigInteger in the specified radix, or {@code radix} is
 *         outside the range from {@link Character#MIN_RADIX} to
 *         {@link Character#MAX_RADIX}, inclusive.
 * @see    Character#digit
 */
BigInt::BigInt(string val, int radix) {
    int cursor = 0;
    len = val.size();

    //to simplify implementation, we do not do input checking here.

    //TODO lianke： still need to translate this java to cpp code.
    // Process first (potentially short) digit group
    int firstGroupLen = len % digitsPerInt[radix];
    if (firstGroupLen == 0)
        firstGroupLen = digitsPerInt[radix];
    string group = val.substring(cursor, cursor += firstGroupLen);
    bytes[BigInt::capacity - 1] = Integer.parseInt(group, radix);

    // Process remaining digit groups
    int superRadix = intRadix[radix];
    int groupVal = 0;
    while (cursor < len) {
        group = val.substring(cursor, cursor += digitsPerInt[radix]);
        groupVal = Integer.parseInt(group, radix);
        destructiveMulAdd(bytes, superRadix, groupVal);
    }
}