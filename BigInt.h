
#include <iostream>
#include <stdexcept>
#include <unistd.h>
#include <cstring>
#include <vector>
#include <string>
#include <bitset>
#define CHECK_BIT(var, pos) (((var) >> (pos)) & 1)

using namespace std;
// TODO we also need to implement BN254, not just for FakeG1 and FakeG2

int uint_msb(uint32_t a)
{
    int n = 4 * sizeof(uint32_t);
    int inf = 0;
    for (int s = n / 2; 0 < s; s /= 2)
    {
        if (a < (uint32_t)1 << n)
        {
            n -= s;
        }
        else
        {
            inf = n;
            n += s;
        }
    }
    return a < (uint32_t)1 << n ? inf : n;
}

int digitsPerInt[] = {0, 0, 30, 19, 15, 13, 11,
                      11, 10, 9, 9, 8, 8, 8, 8, 7, 7, 7, 7, 7, 7, 7, 6, 6, 6, 6,
                      6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 5};

int intRadix[] = {0, 0,
                  0x40000000, 0x4546b3db, 0x40000000, 0x48c27395, 0x159fd800,
                  0x75db9c97, 0x40000000, 0x17179149, 0x3b9aca00, 0xcc6db61,
                  0x19a10000, 0x309f1021, 0x57f6c100, 0xa2f1b6f, 0x10000000,
                  0x18754571, 0x247dbc80, 0x3547667b, 0x4c4b4000, 0x6b5a6e1d,
                  0x6c20a40, 0x8d2d931, 0xb640000, 0xe8d4a51, 0x1269ae40,
                  0x17179149, 0x1cb91000, 0x23744899, 0x2b73a840, 0x34e63b41,
                  0x40000000, 0x4cfa3cc1, 0x5c13d840, 0x6d91b519, 0x39aa400};

class BigInt
{
public:
    uint32_t bytes[16]; // we should have 64 bytes
    static const int capacity = 16;
    static const int num_of_bytes = 64;
    static const int bits_per_word = 32;

    static const long LONG_MASK = 0xffffffffL;

    int len = 0; // number of bytes
    // TODO add modulus here
    BigInt();
    BigInt(string val);
    BigInt(string val, int radix);
    BigInt(uint32_t val);
    BigInt(int val);
    BigInt(const BigInt &val);

    BigInt(unsigned long long val);
    // BigInt(const BigInt& val);

    // Direct assignment
    BigInt &operator=(const BigInt &);

    // Post/Pre - Incrementation
    BigInt &operator++();
    BigInt operator++(int temp);
    BigInt &operator--();
    BigInt operator--(int temp);

    // Addition and Subtraction
    friend BigInt &operator+=(BigInt &, BigInt &);
    friend BigInt operator+(BigInt &, BigInt &);
    friend BigInt operator-(BigInt &, BigInt &);
    friend BigInt &operator-=(BigInt &, const BigInt &);

    // Comparison operators
    friend bool operator==(const BigInt &, const BigInt &);
    friend bool operator!=(const BigInt &, const BigInt &);

    friend bool operator>(const BigInt &, const BigInt &);
    friend bool operator>=(const BigInt &, const BigInt &);
    friend bool operator<(const BigInt &, const BigInt &);
    friend bool operator<=(const BigInt &, const BigInt &);

    // Multiplication and Division
    friend BigInt &operator*=(BigInt &, BigInt &);
    friend BigInt operator*(BigInt &, BigInt &);

    BigInt multiply_with_mod(BigInt &, BigInt &modulus);
    BigInt add_with_mod(BigInt &, BigInt &modulus);
    BigInt minus_with_mod(BigInt &, BigInt &modulus);

    // Modulo
    friend BigInt operator%(BigInt &, BigInt &);
    friend BigInt &operator%=(BigInt &, BigInt &);

    // Power Function
    friend BigInt &operator^=(BigInt &, const BigInt &);
    friend BigInt operator^(BigInt &, const BigInt &);
    friend BigInt pow(BigInt base, int exponent);

    // Read and Write
    void printBinary();
    void printBinaryDebug();

    void printDigits();
    void printAddress();

    void printHex();
    static BigInt ZERO();
    static BigInt ONE();
    int leadingZeros();
    bool isZero();
    bool isOne();
    void shiftLeft(int num_of_bits);
    void shiftRight(int num_of_bits);
    void sbit_helper(int num_of_bits);
    int bitLength();
    bool testBit(int index); // same with the java BigInteger testBit
    BigInt mod(BigInt modulus);
    int getLowestSetBit();
};

int msb(BigInt &a)
{
    int i = 0;
    bool found = 0;
    for (; i < BigInt::capacity; i++)
    {
        if (0 < a.bytes[i])
        {
            found = 1;
            break;
        }
    }
    if (!found)
    {
        return 0;
    }
    return (BigInt::capacity - i - 1) * 8 * sizeof(uint32_t) + uint_msb(a.bytes[i]);
}

BigInt BigInt::ZERO()
{
    return BigInt("0");
}

BigInt BigInt::ONE()
{
    return BigInt("1");
}

BigInt::BigInt()
{
    memset(bytes, 0, num_of_bytes);
}

// TODO lianke implement modulus

BigInt::BigInt(uint32_t val)
{
    memcpy(bytes, &val, num_of_bytes - sizeof(uint32_t));
}

BigInt::BigInt(int val)
{
    // TODO this is wrong
    memcpy(bytes, &val, num_of_bytes - sizeof(int));
}

BigInt::BigInt(unsigned long long val)
{
    memcpy(bytes, &val, num_of_bytes - sizeof(unsigned long long));
}

BigInt::BigInt(const BigInt &val)
{
    memcpy(bytes, &val.bytes, num_of_bytes);
    len = val.len;
}

bool operator>(const BigInt &a, const BigInt &b)
{
    for (int i = 0; i < BigInt::capacity; i++)
    {
        if (a.bytes[i] > b.bytes[i])
        {
            return true;
        }
        else if (a.bytes[i] < b.bytes[i])
        {
            return false;
        }
    }
    return false;
}

void BigInt::shiftLeft(int num_of_bits)
{
    int offset = num_of_bits % BigInt::bits_per_word;
    int num_of_words_shifted = num_of_bits / BigInt::bits_per_word;
    int i = 0;
    if (offset != 0)
    {
        for (; i < BigInt::capacity - num_of_words_shifted - 1; i++)
        {
            bytes[i] = (bytes[i + num_of_words_shifted] << offset) |
                       (bytes[i + num_of_words_shifted + 1] >> (BigInt::bits_per_word - offset));
            // cout << (bytes[i + num_of_words_shifted+1] >> (BigInt::bits_per_word - offset))  << endl;
        }
        // printBinaryDebug();
        bytes[i++] = bytes[i + num_of_words_shifted] << offset;
        for (; i < BigInt::capacity; i++)
        {
            bytes[i] = 0;
        }
    }
    else
    {
        for (; i < BigInt::capacity - num_of_words_shifted - 1; i++)
        {
            bytes[i] = bytes[i + num_of_words_shifted];
        }
        bytes[i++] = bytes[i + num_of_words_shifted];
        for (; i < BigInt::capacity; i++)
        {
            bytes[i] = 0;
        }
    }
}

int LeadingZeros(int x)
{
    const int numIntBits = sizeof(int) * 8; // compile time constant
    // do the smearing
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    // count the ones
    x -= x >> 1 & 0x55555555;
    x = (x >> 2 & 0x33333333) + (x & 0x33333333);
    x = (x >> 4) + x & 0x0f0f0f0f;
    x += x >> 8;
    x += x >> 16;
    return numIntBits - (x & 0x0000003f); // subtract # of 1s from 32
}

int BigInt::bitLength()
{
    int bits = (len - 4) * 8;
    int leading_int = bytes[BigInt::capacity - len / 4];
    return bits + 32 - LeadingZeros(leading_int);
}

int BigInt::leadingZeros()
{
    int res = 0;
    uint32_t zero = 0;
    for (int i = 0; i < BigInt::capacity; i++)
    {
        // cout << "memcmp result=" << memcmp(&bytes[i], &zero, 4)<<endl;
        if (memcmp(&bytes[i], &zero, 4) == 0)
        {
            res += BigInt::bits_per_word;
            // cout <<"continue i=" <<i<<endl;
        }
        else
        {
            res += LeadingZeros(bytes[i]);
            break;
        }
    }
    return res;
}

bool BigInt::isZero()
{
    // TODO lianke test its correctness
    uint32_t testblock[BigInt::capacity];
    memset(testblock, 0, sizeof(testblock));
    return memcmp(testblock, bytes, num_of_bytes) == 0;
}

bool BigInt::isOne()
{
    // TODO lianke test its correctness
    BigInt one("1");
    return *this == one;
}

void BigInt::printBinaryDebug()
{
    uint32_t zero = 0;
    bool leadingZeros = true;
    for (int i = 0; i < capacity; i++)
    {
        if (memcmp(&bytes[i], &zero, sizeof(uint32_t)) == 0 && i != capacity - 1 && leadingZeros)
        {
            continue;
        }
        leadingZeros = false;
        std::bitset<32> tmp(bytes[i]);
        for (int i = BigInt::bits_per_word - 1; i >= 0; i--)
        {
            cout << tmp[i];
            if (i % 8 == 0)
            {
                cout << "|";
            }
        }
    }
    printf("\n");
    return;
}

void BigInt::printBinary()
{
    uint32_t zero = 0;
    bool leadingZeros = true;
    for (int i = 0; i < capacity; i++)
    {
        if (memcmp(&bytes[i], &zero, sizeof(uint32_t)) == 0 && i != capacity - 1 && leadingZeros)
        {
            continue;
        }
        leadingZeros = false;
        std::bitset<32> tmp(bytes[i]);
        cout << tmp;
    }
    printf("\n");
    return;
}

void BigInt::printHex()
{
    for (int i = 0; i < capacity; i++)
    {
        printf("%X", bytes[i]);
    }
    printf("\n");
    return;
}

BigInt::BigInt(string val) : BigInt(val, 10)
{
}

void BigInt::printAddress()
{
    printf("0%x\n", bytes);
    return;
}

bool BigInt::testBit(int n)
{
    int byte_index = BigInt::capacity - 1 - (n / BigInt::bits_per_word);
    int byte_offset = n % BigInt::bits_per_word;
    // printf("%hhx %hhx ", bytes[byte_index], byte_offset);

    return CHECK_BIT(bytes[byte_index], byte_offset);
    // return (bytes[n >> 5] & (1 << (n & 31))) != 0;
}

BigInt &BigInt::operator=(const BigInt &a)
{
    memcpy(bytes, a.bytes, BigInt::num_of_bytes);
    return *this;
}

bool operator==(const BigInt &a, const BigInt &b)
{
    return memcmp(a.bytes, b.bytes, BigInt::num_of_bytes) == 0;
}

bool operator<(const BigInt &a, const BigInt &b)
{
    for (int i = 0; i < BigInt::capacity; i++)
    {
        if (a.bytes[i] == 0 && b.bytes[i] == 0)
        {
            continue;
        }
        return a.bytes[i] < b.bytes[i];
    }
}

BigInt &operator%=(BigInt &a, BigInt &b)
{
    a = a % b;
}

BigInt operator%(BigInt &a, BigInt &b)
{
    BigInt div_result;
    BigInt mod_result;
    int msb_a = msb(a);
    int msb_b = msb(b);
    if (msb_a < msb_b)
    {
        mod_result = a;
    }
    else
    {
        int xbit = msb_a - msb_b;
        // cout << "msb diff = " << xbit << endl;
        BigInt xobj = a;
        BigInt xdiv = b;
        xdiv.shiftLeft(xbit);
        for (int i = 0; i <= xbit; ++i)
        {
            // cout << "xobj binary = ";
            // xobj.printBinary();
            // cout << "xdiv binary = ";
            // xdiv.printBinary();
            // cout << endl;
            if (xobj > xdiv)
            {
                div_result.sbit_helper(xbit - i);
                xobj = xobj - xdiv;
            }
            xdiv.shiftRight(1);
        }
        mod_result = xobj;
    }
    return mod_result;
}

uint32_t uint_sub_helper(uint32_t a, uint32_t b, bool &carry)
{
    uint32_t cx = carry ? 1 : 0;
    uint32_t c = a - b - cx;
    carry = carry ? a <= c : a < c;
    return c;
}

BigInt operator-(BigInt &a, BigInt &b)
{
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;

    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        result.bytes[i] = uint_sub_helper(a.bytes[i], b.bytes[i], carry);
    }
    if (carry)
    {
        printf("-------------------BigInt subtraction failure detected!!----------------\n");
    }
    return result;
}

BigInt BigInt::minus_with_mod(BigInt &b, BigInt &modulus)
{
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;
    BigInt a = *this;
    if (a < b)
    {
        a = a + modulus;
    }

    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        result.bytes[i] = uint_sub_helper(a.bytes[i], b.bytes[i], carry);
    }

    a = a % modulus;
    return result;
}

BigInt operator+(BigInt &a, BigInt &b)
{
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;
    // lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        // cout << a.bytes[i] << " " << b.bytes[i] << " " <<temp << endl;
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        result.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigInt::bits_per_word != 0);
    }
    return result;
}

BigInt BigInt::add_with_mod(BigInt &b, BigInt &modulus)
{
    BigInt result;
    uint64_t temp = 0;
    bool carry = false;
    // lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        // cout << a.bytes[i] << " " << b.bytes[i] << " " <<temp << endl;
        temp = (uint64_t)bytes[i] + b.bytes[i] + carry;
        result.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigInt::bits_per_word != 0);
    }
    result = result % modulus;
    return result;
}

BigInt &operator+=(BigInt &a, const BigInt &b)
{
    uint64_t temp = 0;
    bool carry = false;
    // lianke: we only use the lower half. the higher half is for storing larger multiplication results.
    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        temp = (uint64_t)a.bytes[i] + b.bytes[i] + carry;
        a.bytes[i] = (uint32_t)temp;
        carry = (temp >> BigInt::bits_per_word != 0);
    }
    return a;
}

BigInt FqModulusParameter("1532495540865888858358347027150309183618765510462668801");

/* Returns this^exponent with long exponent */
// TODO lianke : check its correctness
BigInt pow(BigInt base, int exponent, BigInt modulus)
{
    BigInt value = base;
    BigInt result("1");
    int currentExponent = exponent;
    while (currentExponent > 0)
    {
        // cout << "currentExponent=" << currentExponent <<endl;
        //  result.printBinary();
        //  value.printBinary();
        if (currentExponent % 2 == 1)
        {
            result = result * value;
            result %= modulus;
        }
        value = value * value;
        value %= modulus;
        currentExponent >>= 1;
    }
    return result;
}

uint32_t uint_add_helper(uint32_t a, uint32_t b, bool &carry)
{
    uint64_t result = (uint64_t)a + (uint64_t)b;
    carry = (result >> 32) != 0;
    return result;
}

BigInt BigInt::multiply_with_mod(BigInt &b, BigInt &modulus)
{
    BigInt result;
    BigInt carry;
    // cout << "start of a multiplication" << endl;

    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        for (int j = BigInt::capacity - 1; j >= 0; j--)
        {
            int k = i + j + 1 - BigInt::capacity;
            if (k >= 0)
            {
                uint64_t kmul = (uint64_t)bytes[i] * b.bytes[j];
                uint32_t kmul_lower = (uint32_t)kmul;
                uint32_t kmul_higher = kmul >> 32;
                bool carry1 = false;
                bool carry2 = false;
                result.bytes[k] = uint_add_helper(result.bytes[k], kmul_lower, carry1);
                if (k >= 1)
                {
                    if (carry1)
                    {
                        carry.bytes[k - 1]++;
                    }
                    result.bytes[k - 1] = uint_add_helper(result.bytes[k - 1], kmul_higher, carry2);
                    if (carry2 && k >= 2)
                    {
                        carry.bytes[k - 2]++;
                    }
                }
            }
        }
    }

    BigInt final_res = result + carry;
    final_res = final_res % modulus;
    return final_res;
}

BigInt &operator*=(BigInt &a, BigInt &b)
{

    a = a * b;
    return a;
}

BigInt operator*(BigInt &a, BigInt &b)
{

    BigInt result;
    BigInt carry;
    // cout << "start of a multiplication" << endl;

    for (int i = BigInt::capacity - 1; i >= 0; i--)
    {
        for (int j = BigInt::capacity - 1; j >= 0; j--)
        {
            int k = i + j + 1 - BigInt::capacity;
            if (k >= 0)
            {
                uint64_t kmul = (uint64_t)a.bytes[i] * b.bytes[j];
                uint32_t kmul_lower = (uint32_t)kmul;
                uint32_t kmul_higher = kmul >> 32;
                bool carry1 = false;
                bool carry2 = false;
                result.bytes[k] = uint_add_helper(result.bytes[k], kmul_lower, carry1);
                if (k >= 1)
                {
                    if (carry1)
                    {
                        carry.bytes[k - 1]++;
                    }
                    result.bytes[k - 1] = uint_add_helper(result.bytes[k - 1], kmul_higher, carry2);
                    if (carry2 && k >= 2)
                    {
                        carry.bytes[k - 2]++;
                    }
                }
            }
        }
    }

    BigInt final_res = result + carry;

    return final_res;
}

void BigInt::shiftRight(int num_of_bits)
{
    // TODO lianke this shift is wrong
    int offset = num_of_bits % BigInt::bits_per_word;
    int num_of_words_shifted = num_of_bits / BigInt::bits_per_word;
    int i = BigInt::capacity - 1;
    if (offset != 0)
    {
        for (; i >= num_of_words_shifted + 1; i--)
        {
            bytes[i] = (bytes[i - num_of_words_shifted] >> offset) |
                       (bytes[i - num_of_words_shifted - 1] << (BigInt::bits_per_word - offset));
        }
        bytes[i--] = bytes[i - num_of_words_shifted] >> offset;
        for (; i >= 0; i--)
        {
            bytes[i] = 0;
        }
    }
    else
    {
        for (; i >= num_of_words_shifted + 1; i--)
        {
            bytes[i] = bytes[i - num_of_words_shifted];
        }
        bytes[i--] = bytes[i - num_of_words_shifted];
        for (; i >= 0; i--)
        {
            bytes[i] = 0;
        }
    }
}

void BigInt::sbit_helper(int bit)
{
    int dat_byte = bit / BigInt::bits_per_word;
    int dat_bit = bit % BigInt::bits_per_word;
    bytes[dat_byte] |= 1 << dat_bit;
}

BigInt mod(BigInt &a, BigInt &b)
{
    BigInt div_result;
    BigInt mod_result;
    int msb_a = msb(a);
    int msb_b = msb(b);
    if (msb_a < msb_b)
    {
        mod_result = a;
    }
    else
    {
        int xbit = msb_a - msb_b;
        // cout << "msb diff = " << xbit << endl;
        BigInt xobj = a;
        BigInt xdiv = b;
        xdiv.shiftLeft(xbit);
        for (int i = 0; i <= xbit; ++i)
        {
            // cout << "xobj binary = ";
            // xobj.printBinary();
            // cout << "xdiv binary = ";
            // xdiv.printBinary();
            // cout << endl;
            if (xobj > xdiv)
            {
                div_result.sbit_helper(xbit - i);
                xobj = xobj - xdiv;
            }
            xdiv.shiftRight(1);
        }
        mod_result = xobj;
    }
    return mod_result;
}

// Multiply x array times word y in place, and add word z
void destructiveMulAdd(uint32_t x[], int y, int z)
{
    // Perform the multiplication word by word
    long ylong = y & BigInt::LONG_MASK;
    long zlong = z & BigInt::LONG_MASK;
    int len = BigInt::capacity;

    unsigned long long product = 0;
    unsigned long carry = 0;
    for (int i = len - 1; i >= 0; i--)
    {
        product = ylong * (x[i] & BigInt::LONG_MASK) + carry;
        x[i] = (uint32_t)product;
        carry = product >> 32;
    }

    // Perform the addition
    unsigned long sum = (x[len - 1] & BigInt::LONG_MASK) + zlong;
    x[len - 1] = (uint32_t)sum;
    carry = sum >> 32;
    for (int i = len - 2; i >= 0; i--)
    {
        sum = (x[i] & BigInt::LONG_MASK) + carry;
        x[i] = (uint32_t)sum;
        carry = sum >> 32;
    }
}

BigInt::BigInt(string val, int radix)
{
    int cursor = 0;
    len = val.size();
    memset(bytes, 0, BigInt::num_of_bytes);
    // to simplify implementation, we do not do input checking here.

    // Process first (potentially short) digit group
    int firstGroupLen = len % digitsPerInt[radix];
    if (firstGroupLen == 0)
        firstGroupLen = digitsPerInt[radix];
    string group = val.substr(cursor, firstGroupLen);
    cursor += firstGroupLen;
    bytes[BigInt::capacity - 1] = stoi(group);
    // cout << group << endl;

    // Process remaining digit groups
    int superRadix = intRadix[radix];
    int groupVal = 0;
    while (cursor < len)
    {
        // cout << cursor << " " << cursor + digitsPerInt[radix] << endl;
        group = val.substr(cursor, digitsPerInt[radix]);
        cursor += digitsPerInt[radix];
        // cout << "read in :" << group << endl;
        groupVal = stoi(group);
        destructiveMulAdd(bytes, superRadix, groupVal);
    }
}
