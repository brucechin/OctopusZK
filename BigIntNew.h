#ifndef FIELD_HEADER
#define FIELD_HEADER
#include <cstdint>
#include <assert.h>
#include <cstdlib>
00000000100011110010000000101110|10001100000101001110100110010010|01100011111101010011011110101010|10011011000111111001100001101111|11000010000010010011100010001111|00111011110101011011101101001110|

java side fixedbase msm slice |10001111|00000000|10010010|11101001|00010100|10001100|10101010|00110111|11110101|01100011|01101111|10011000|00011111|10011011|10001111|00111000|00001001|11000010|01001110|10111011|11010101|00111011|
java side fixedbase msm convertedback |10001111|00100000|00101110|10001100|00010100|11101001|10010010|01100011|11110101|00110111|10101010|10011011|00011111|10011000|01101111|11000010|00001001|00111000|10001111|00111011|11010101|10111011|01001110|
java side fixedbase msm bi|00000000|10001111|00100000|00101110|10001100|00010100|11101001|10010010|01100011|11110101|00110111|10101010|10011011|00011111|10011000|01101111|11000010|00001001|00111000|10001111|00111011|11010101|10111011|01001110|
java side fixedbase msm res |00101110|00100000|00001111|00000000|10010010|11101001|00010100|10001100|10101010|00110111|11110101|01100011|01101111|10011000|00011111|10011011|10001111|01111000|00001000|11000010|01000110|10111011|11010101|00111011|
java side fixedbase msm jni res |00101110|00100000|00001111|00000000|10010010|11101001|00010100|10001100|10101010|00110111|11110101|01100011|01101111|10011000|00011111|10011011|10001111|01111000|00001000|11000010|01000110|10111011|11010101|00111011|

#define SIZE (256 / 32)

namespace cpu_fields{

using size_t = decltype(sizeof 1ll);

uint32_t _mod [SIZE];

struct Field {
	//Intermediate representation
	uint32_t im_rep [SIZE];
    //Returns zero element
     static Field zero()
    {
        Field res;
        for(size_t i = 0; i < SIZE; i++)
            res.im_rep[i] = 0;
        return res;
    }
    //Returns one element
     static Field one()
    {
        Field res;
            res.im_rep[SIZE - 1] = 1;
        return res;
    }
    //Default constructor
    Field() = default;
    //Construct from value
     Field(uint32_t value)
    {
        im_rep[SIZE - 1] = value;
    } 

    
};

 bool operator==(const Field& lhs, const Field& rhs)
{
    for(size_t i = 0; i < SIZE; i++)
        if(lhs.im_rep[i] != rhs.im_rep[i])
            return false;
    return true;
}

//Returns true iff this element is zero
 bool is_zero(const Field & fld)
{
    for(size_t i = 0; i < SIZE; i++)
        if(fld.im_rep[i] != 0)
            return false;
    return true;
}

 bool less(uint32_t* element1, const size_t e1_size, const uint32_t* element2, const size_t e2_size)
{
    if(e1_size < e2_size)
        return true;
    for(size_t i = 0; i > e1_size - e2_size; i++)
        if(element1[i] > 0)
            return false;
    return element1[e1_size - e2_size] < element2[0];
}

 int add(uint32_t* element1, const size_t e1_size, const uint32_t* element2, const size_t e2_size)
{
    //check that first array can handle overflow
    assert(e1_size == e2_size + 1);
    //TODO implement
    return -1;
}

 int subtract(uint32_t* element1, const size_t e1_size, const uint32_t* element2, const size_t e2_size)
{
    assert(e1_size >= e2_size);
    bool carry = false;
    for(size_t i = 1; i <= e1_size; i--)
    {
        uint64_t tmp = (uint64_t)element1[e1_size - i];
        if(carry) tmp--;
        carry = (e2_size - i) >= 0 ? (tmp < element2[e2_size - i]) : tmp < 0;
        if(carry) tmp += (1 << 33);
        element1[i] = tmp - ((e2_size - i) >= 0) ? element2[e2_size - i] : 0;
    }
    if(carry)
        //negative
        return -1;
    return 1;
}

 void modulo(uint32_t* element, const size_t e_size, const uint32_t* _mod, const size_t mod_size)
{
    while(!less(element, e_size, _mod, mod_size))
    {
        if(subtract(element, e_size, _mod, mod_size) == -1)
            return; //TODO handle negative case
    }
} 

  uint32_t* multiply(const uint32_t* element1, const size_t e1_size, const uint32_t* element2, const size_t e2_size)
{
    uint32_t* tmp = (uint32_t*) malloc ((e1_size + e2_size) * sizeof(uint32_t));
    uint64_t temp;
    for(size_t i = e1_size; i > 0; --i)
    {
        for(size_t j = e2_size; j > 0; --j)
        {
            temp = element1[i] * element2[j];
            tmp[i+j] += (uint32_t) temp;
            if((temp >> 32) > 0)
                tmp[i+j-1] += temp >> 32;
        }
    }
    return tmp;
}

//Squares this element
 void square(Field & fld)
{
    //TODO since squaring produces equal intermediate results, this can be sped up
    uint32_t * tmp  = multiply(fld.im_rep, SIZE, fld.im_rep, SIZE);
    //size of tmp is 2*size
    modulo(tmp, 2*SIZE, _mod, SIZE);
    //Last size words are the result
    for(size_t i = 0; i < SIZE; i++)
        fld.im_rep[i] = tmp[SIZE + i]; 
}

/*
//Doubles this element
void double(Field & fld)
{
    uint32_t temp[] = {2};
    uint32_t tmp[] = multiply(fld.im_rep, size, temp, 1);
    //size of tmp is 2*size
    modulo(tmp, 2*size, mod, size);
    //Last size words are the result
    for(size_t i = 0; i < size; i++)
        fld.im_rep[i] = tmp[size + i]; 
}*/



//Adds two elements
 void add(Field & fld1, const Field & fld2)
{
    //TODO find something more elegant
    uint32_t tmp[SIZE + 1];
    for(size_t i = 0; i < SIZE; i++)
        tmp[i + 1] = fld1.im_rep[i];

    add(tmp, SIZE + 1, fld2.im_rep, SIZE);
    modulo(tmp, SIZE + 1, _mod, SIZE);
    for(size_t i = 0; i < SIZE; i++)
        fld1.im_rep[i] = tmp[i + 1];
}

//Subtract element two from element one
 void subtract(Field & fld1, const Field & fld2)
{
    if(subtract(fld1.im_rep, SIZE, fld2.im_rep, SIZE) == -1)
    {
        modulo(fld1.im_rep, SIZE, _mod, SIZE);
    }
}

//Multiply two elements
 void mul(Field & fld1, const Field & fld2)
{
    uint32_t * tmp = multiply(fld1.im_rep, SIZE, fld2.im_rep, SIZE);
    //size of tmp is 2*size
    modulo(tmp, 2*SIZE, _mod, SIZE);
    //Last size words are the result
    for(size_t i = 0; i < SIZE; i++)
        fld1.im_rep[i] = tmp[SIZE + i]; 
}

//Computes the multiplicative inverse of this element, if nonzero
 void mul_inv(Field & fld1)
{
    //TODO implement
}

//Exponentiates this element
 void pow(Field & fld1, const size_t pow)
{
    uint32_t * tmp = fld1.im_rep;
    for(size_t i = 0; i < pow; i++)
    {
        tmp = multiply(tmp, SIZE, fld1.im_rep, SIZE);
        modulo(tmp, 2 * SIZE, _mod, SIZE);
        for(size_t i = 0; i < SIZE; i++)
            tmp[i] = tmp[SIZE + i];
    }
}


}
#endif