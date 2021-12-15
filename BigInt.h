
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
    BigInt(const BigInt& val);


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
    friend BigInt operator-( BigInt &,  BigInt &);
    friend BigInt &operator-=(BigInt &, const BigInt &);
 
    //Comparison operators
    friend bool operator==(const BigInt &, const BigInt &);
    friend bool operator!=(const BigInt &, const BigInt &);
 
    friend bool operator>(const BigInt &, const BigInt &);
    friend bool operator>=(const BigInt &, const BigInt &);
    friend bool operator<(const BigInt &, const BigInt &);
    friend bool operator<=(const BigInt &, const BigInt &);

    //Multiplication and Division
    friend BigInt &operator*=(BigInt &,  BigInt &);
    friend BigInt operator*( BigInt &,  BigInt &);

    //Modulo
    friend BigInt operator%( BigInt &,  BigInt &);
    friend BigInt &operator%=(BigInt &,  BigInt &);

    //Power Function
    friend BigInt &operator^=(BigInt &,const BigInt &);
    friend BigInt operator^(BigInt &, const BigInt &);
    friend BigInt pow(BigInt base, int exponent);

    //Read and Write
    void printBinary();
    void printDigits();
    void printAddress();

    void printHex();
    static BigInt ZERO();
    static BigInt ONE();
    int leadingZeros();
    bool isZero();
    bool isOne();
    int bitLength();
    bool testBit(int index); //same with the java BigInteger testBit
    BigInt mod(BigInt modulus);
    int getLowestSetBit() ;
};




BigInt BigInt::ZERO(){
    return BigInt("0");
}

BigInt BigInt::ONE(){
    return BigInt("1");
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
    //TODO this is wrong
    memcpy(bytes, &val, num_of_bytes - sizeof(int));
}

BigInt::BigInt(unsigned long long val){
    memcpy(bytes, &val, num_of_bytes - sizeof(unsigned long long ));
}

BigInt::BigInt(const BigInt& val){
    memcpy(bytes, &val.bytes, num_of_bytes);
    len  = val.len;
}








int LeadingZeros(int x)
{
    const int numIntBits = sizeof(int) * 8; //compile time constant
    //do the smearing
    x |= x >> 1; 
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    //count the ones
    x -= x >> 1 & 0x55555555;
    x = (x >> 2 & 0x33333333) + (x & 0x33333333);
    x = (x >> 4) + x & 0x0f0f0f0f;
    x += x >> 8;
    x += x >> 16;
    return numIntBits - (x & 0x0000003f); //subtract # of 1s from 32
}

int BigInt::bitLength(){
    int bits = (len - 4) * 8;
    int leading_int = bytes[BigInt::capacity - len/4];
    return bits + 32 - LeadingZeros(leading_int);
}

int BigInt::leadingZeros(){
    int res = 0;
    uint32_t zero = 0;
    for(int i = 0; i < BigInt::capacity; i++){
        //cout << "memcmp result=" << memcmp(&bytes[i], &zero, 4)<<endl;
        if(memcmp(&bytes[i], &zero, 4) == 0){
            res += BigInt::bits_per_word;
            //cout <<"continue i=" <<i<<endl; 
        }else{
            res +=LeadingZeros(bytes[i]);
            break;
        }
    }
    return res;

}

bool BigInt::isZero(){
    //TODO lianke test its correctness
    uint32_t testblock[BigInt::capacity];
    memset(testblock, 0, sizeof(testblock));
    return memcmp(testblock, bytes, num_of_bytes) == 0;
}

bool BigInt::isOne(){
    //TODO lianke test its correctness
    BigInt one("1");
    return *this == one;
}



void BigInt::printBinary(){
    uint32_t zero = 0;
    bool leadingZeros = true;
    for (int i = 0; i < capacity; i++){
        if(memcmp(&bytes[i], &zero, sizeof(uint32_t)) == 0 && i !=capacity-1 && leadingZeros){ continue; }
        leadingZeros = false; 
        std::bitset<32> tmp(bytes[i]);
        cout << tmp << "|";
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
    int byte_index = BigInt::capacity - 1 - (n / BigInt::bits_per_word); 
    int byte_offset =  n % BigInt::bits_per_word;
    //printf("%hhx %hhx ", bytes[byte_index], byte_offset);
    
    return CHECK_BIT(bytes[byte_index], byte_offset);
    //return (bytes[n >> 5] & (1 << (n & 31))) != 0;
}

BigInt &BigInt::operator=(const BigInt &a){
    memcpy(bytes, a.bytes, BigInt::num_of_bytes);
    return *this;
}

bool operator==(const BigInt &a, const BigInt &b){
    return memcmp(a.bytes, b.bytes, BigInt::num_of_bytes) == 0;
}

bool operator<(const BigInt &a, const BigInt &b){
    for(int i = 0; i < BigInt::capacity; i++){
        if(a.bytes[i] == 0 && b.bytes[i] == 0){ 
            continue;
        }
        return a.bytes[i] < b.bytes[i];
    }
}

BigInt &operator%=(BigInt &a,  BigInt &b){
    //TODO lianke optimize this func
    int count = 0;
    while(!(a < b)){
        //cout << "mod index=" << count++ <<endl;
        //a.printBinary();
        //b.printBinary();
        int leadingZerosA = a.leadingZeros();
        int leadingZerosB = b.leadingZeros();
        int leadingZeroDiff = leadingZerosB - leadingZerosA;
        cout << leadingZerosA << " " <<leadingZerosB << " " << leadingZeroDiff <<endl;
        if(leadingZeroDiff > 0){ 
            BigInt multiplier("0");
            int index = BigInt::capacity - 1 - (leadingZeroDiff / BigInt::bits_per_word);
            uint32_t offset = 1 << ((leadingZeroDiff - 1) % BigInt::bits_per_word);
            //cout << "index=" << index << " offset=" <<offset <<endl;
            //memcpy(&multiplier.bytes[index], &offset, sizeof(uint32_t));
            multiplier.bytes[index] = offset;
            multiplier.printBinary();
            BigInt tmp = b;
            cout <<"previsous tmp" <<endl;
            tmp.printBinary();
            tmp = tmp * multiplier;
            cout << "multiplied tmp:" << endl;
            tmp.printBinary();

            a = a - tmp;
        }else{
            a = a - b;
        }

    }
    return a;
}

BigInt operator%( BigInt &a,  BigInt &b){
    BigInt temp;
    temp = a;
    temp %= b;
    return temp;
}

BigInt operator-(BigInt &a,  BigInt& b) {
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;
    //lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for(int i = BigInt::capacity - 1; i >= 0; i--) {
        //cout << a.bytes[i] << " " << b.bytes[i] << " " <<temp << endl;
        temp = (uint64_t)a.bytes[i];
        if(carry){
            temp--;
        }
        carry = temp < (uint64_t)b.bytes[i];
        if(carry){
            temp += (1 << 32);
        }
        result.bytes[i] = (uint32_t)(temp - (uint64_t)b.bytes[i]);
    }
    if(carry){
        printf("-------------------BigInt subtraction failure detected!!----------------\n");
    }
    return result;
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
    for(int i = BigInt::capacity - 1; i >= 0; i--) {
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        a.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigInt::bits_per_word != 0);
    }
    return a;
}


BigInt FqModulusParameter("1532495540865888858358347027150309183618765510462668801");

/* Returns this^exponent with long exponent */
BigInt pow(BigInt base, int exponent) {
    BigInt value = base;
    BigInt result("1");
    int currentExponent = exponent;
    while (currentExponent > 0) {
        cout << "currentExponent=" << currentExponent <<endl;
        result = result * value;

        result.printBinary();
        value.printBinary();
        result %= FqModulusParameter;

        // if (currentExponent % 2 == 1) {
        //     result = result * value;

        //     result %= FqModulusParameter;
        // }
        value  = value * value;
        value %= FqModulusParameter;
        currentExponent >>= 1;
    }
    return result;
}


//     /**
//      * Compare the magnitude of two MutableBigIntegers. Returns -1, 0 or 1
//      * as this MutableBigInteger is numerically less than, equal to, or
//      * greater than <tt>b</tt>.
//      */
//     int compare(BigInt b) {
//         int blen = b.len;
//         if (len < blen) 
//             return -1;
//         if (len > blen)
//            return 1;

//         // Add Integer.MIN_VALUE to make the comparison act as unsigned integer
//         // comparison.
//         int bval[BigInt::capacity] = b.bytes;
//         for (int i = BigInt::capacity - 1; i >=0; i--) {
//             int b1 = value[i] + 0x80000000;
//             int b2 = bval[i]  + 0x80000000;
//             if (b1 < b2)
//                 return -1;
//             if (b1 > b2)
//                 return 1;
//         }
//         return 0;
//     }


//     int numberOfTrailingZeros(int i) {
//         // HD, Figure 5-14
//         int y;
//         if (i == 0) return 32;
//         int n = 31;
//         y = i <<16; if (y != 0) { n = n -16; i = y; }
//         y = i << 8; if (y != 0) { n = n - 8; i = y; }
//         y = i << 4; if (y != 0) { n = n - 4; i = y; }
//         y = i << 2; if (y != 0) { n = n - 2; i = y; }
//         return n - ((i << 1) >>> 31);
//     }

//    /**
//      * Return the index of the lowest set bit in this MutableBigInteger. If the
//      * magnitude of this MutableBigInteger is zero, -1 is returned.
//      */
//     int BigInt::getLowestSetBit() {
//         if (len == 0)
//             return -1;
//         int j, b;
//         for (j=BitInt::capacity-1; (j >= 0) && (bytes[j] == 0); j--)
//             ;
//         b = value[j];
//         if (b == 0)
//             return -1;
//         return ((BigInt::capacity-1-j)<<5) + numberOfTrailingZeros(b);
//     }


//     BigInt divideKnuth(BigInt b, BigInt quotient, bool needRemainder) {

//         // Dividend is zero
//         if (intLen == 0) {
//             quotient.intLen = quotient.offset = 0;
//             return needRemainder ? new MutableBigInteger() : null;
//         }

//         int cmp = compare(b);
//         // Dividend less than divisor
//         if (cmp < 0) {
//             quotient.intLen = quotient.offset = 0;
//             return needRemainder ? BigInt(this): null;
//         }


//         quotient.clear();
//         // Special case one word divisor
//         if (b.intLen == 1) {
//             int r = divideOneWord(b.value[b.offset], quotient);
//             if(needRemainder) {
//                 if (r == 0)
//                     return new MutableBigInteger();
//                 return new MutableBigInteger(r);
//             } else {
//                 return null;
//             }
//         }

//         // Cancel common powers of two if we're above the KNUTH_POW2_* thresholds
//         if (intLen >= KNUTH_POW2_THRESH_LEN) {
//             int trailingZeroBits = Math.min(getLowestSetBit(), b.getLowestSetBit());
//             if (trailingZeroBits >= KNUTH_POW2_THRESH_ZEROS*32) {
//                 MutableBigInteger a = new MutableBigInteger(this);
//                 b = new MutableBigInteger(b);
//                 a.rightShift(trailingZeroBits);
//                 b.rightShift(trailingZeroBits);
//                 MutableBigInteger r = a.divideKnuth(b, quotient);
//                 r.leftShift(trailingZeroBits);
//                 return r;
//             }
//         }

//         return divideMagnitude(b, quotient, needRemainder);
//     }



BigInt &operator*=(BigInt &a, BigInt& b){

    //TODO lianke implement Karatsuba multiplication.
    uint64_t temp[BigInt::capacity * 2] = {0}; 
    for(int i = BigInt::capacity - 1; i >= 0; i--){
        for(int j = BigInt::capacity - 1; j>=0; j--){
            temp[i + j  + 1] += (uint64_t)a.bytes[i] * b.bytes[j];
        }
    }
    uint64_t tmp = 0;
    uint64_t carry = 0;
    for(int i = BigInt::capacity - 1; i >= 0; i--){
        tmp = temp[BigInt::capacity + i] + carry;
        a.bytes[i] = (uint32_t)tmp;
        carry = (tmp >> BigInt::bits_per_word);
    }

    a = a % FqModulusParameter;

    return a;
}


BigInt operator*(BigInt &a, BigInt& b){

    BigInt result;
    //TODO lianke implement Karatsuba multiplication.
    uint64_t temp[BigInt::capacity * 2] = {0}; 
    for(int i = BigInt::capacity - 1; i >= 0; i--){
        for(int j = BigInt::capacity - 1; j>= 0; j--){
            temp[i + j  + 1] += (uint64_t)a.bytes[i] * b.bytes[j];
        }
    }

    uint64_t tmp = 0;
    uint64_t carry = 0;
    for(int i = BigInt::capacity - 1; i >= 0; i--){
        tmp = temp[BigInt::capacity + i] + carry;
        result.bytes[i] = (uint32_t)tmp;
        carry = (tmp >> BigInt::bits_per_word);
    }
    result = result % FqModulusParameter;

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



