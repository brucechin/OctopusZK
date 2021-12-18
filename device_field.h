#pragma once
#include <cstdint>

#ifndef DEBUG
#include <cuda.h>
#include <cuda_runtime.h>

#define cu_fun __host__ __device__
#else

#define cu_fun
#include <cstdio>
#include <cassert>

#endif

#define SIZE (512 / 32)
#if defined(__CUDA_ARCH__)
#define _clz __clz
#else
#define _clz __builtin_clz
#endif


#define CHECK_BIT(var, pos) ((var) & (1 << (pos)))

#define m_inv 4294967295L

using size_t = decltype(sizeof 1ll);

#ifndef DEBUG
__constant__
#endif


    //_mod for inverted representation
    //TODO Lianke we need to check its endian correctness
// const uint32_t _mod[SIZE] = {0,0,0,0, 
//                             0,0,0,0,
//                             0, 0,1048576,0,
//                             0,0,6144,1};
const uint32_t _mod[SIZE] = {1,6144,0,0, 
                            0,1048576,0,0,
                            0, 0,0,0,
                            0,0,0,0};
struct Scalar
{
    cu_fun void add(Scalar &fld1, const Scalar &fld2) const;
    cu_fun void mul(Scalar &fld1, const Scalar &fld2, const uint32_t *mod) const;
    cu_fun void subtract(Scalar &fld1, const Scalar &fld2) const;
    cu_fun void pow(Scalar &fld1, const uint32_t pow) const;
    static const int capacity = 16; 
    static const int bits_per_word = 32; 
    static const int num_of_bytes = 64;
    //Intermediate representation
    uint32_t im_rep[SIZE] = {0};
    //Returns zero element
    cu_fun static Scalar zero()
    {
        Scalar res;
        for (size_t i = 0; i < SIZE; i++)
            res.im_rep[i] = 0;
        return res;
    }

    //Returns one element
    cu_fun static Scalar one()
    {
        Scalar res;
        res.im_rep[0] = 1;
        return res;
    }
    //Default constructor
    cu_fun Scalar() = default;
    //Construct from value
    cu_fun Scalar(const uint32_t value)
    {
        im_rep[0] = value;
    }

    cu_fun Scalar(const uint32_t *value)
    {
        for (size_t i = 0; i < SIZE; i++)
            im_rep[i] = value[i];
    }

    //Returns true iff this element is zero
    cu_fun bool is_zero() const
    {
        for (size_t i = 0; i < SIZE; i++)
            if (this->im_rep[i] != 0)
                return false;
        return true;
    }

    cu_fun Scalar operator*(const Scalar &rhs) const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        mul(s, rhs, _mod);
        return s;
    }

    cu_fun Scalar operator+(const Scalar &rhs) const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        add(s, rhs);
        return s;
    }

    cu_fun Scalar operator-(const Scalar &rhs) const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        subtract(s, rhs);
        return s;
    }

    cu_fun Scalar operator-() const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        subtract(s, *this);
        return s;
    }

    cu_fun Scalar operator^(const uint32_t &rhs) const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        //TODO lianke fix this pow
        pow(s, rhs);
        return s;
    }

    cu_fun bool operator==(const Scalar &rhs) const
    {
        for (size_t i = 0; i < SIZE; i++)
            if (rhs.im_rep[i] != this->im_rep[i])
                return false;
        return true;
    }
    /*
    cu_fun Scalar operator=(const Scalar &rhs) const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = rhs.im_rep[i];
        return s;
    }*/

    cu_fun Scalar square() const
    {
        Scalar s;
        for (size_t i = 0; i < SIZE; i++)
            s.im_rep[i] = this->im_rep[i];
        mul(s, *this, _mod);
        return s;
    }

    cu_fun static Scalar shuffle_down(unsigned mask, Scalar val, unsigned offset)
    {
        Scalar result;
        for (size_t i = 0; i < SIZE; i++)
#if defined(__CUDA_ARCH__)
            result.im_rep[i] = __shfl_down_sync(mask, val.im_rep[i], offset);
#else
            result.im_rep[i] = val.im_rep[i];
#endif
        return result;
    }

    cu_fun static void print(Scalar f)
    {
        for (size_t i = 0; i < SIZE; i++)
            printf("%u, ", f.im_rep[i]);
        printf("\n");
    }

    static void testEquality(Scalar f1, Scalar f2)
    {
        for (size_t i = 0; i < SIZE; i++)
            if (f1.im_rep[i] != f2.im_rep[i])
            {
                printf("Missmatch: \n");
                print(f1);
                print(f2);
                assert(!"missmatch");
            }
    }
};

cu_fun long idxOfLNZ(const Scalar &fld);
cu_fun bool hasBitAt(const Scalar &fld, long index);



using size_t = decltype(sizeof 1ll);

cu_fun bool operator==(const Scalar &lhs, const Scalar &rhs)
{
    for (size_t i = 0; i < SIZE; i++)
        if (lhs.im_rep[i] != rhs.im_rep[i])
            return false;
    return true;
}

cu_fun uint32_t clz(const uint32_t *element, const size_t e_size)
{
    uint32_t lz = 0;
    uint32_t tmp;
    for (size_t i = e_size; i > 0; i--)
    {
        if (element[i] == 0)
            tmp = 32;
        else
            tmp = _clz(element[i]);
        lz += tmp;
        if (tmp < 32)
            break;
    }
    return lz;
}

cu_fun long idxOfLNZ(const Scalar &fld)
{
    return SIZE - clz(fld.im_rep, SIZE);
}

cu_fun bool hasBitAt(const Scalar &fld, long index)
{
    long idx1 = index % 32;
    long idx2 = index / 32;
    return CHECK_BIT(fld.im_rep[idx2], idx1) != 0;
}

//Returns true if the first element is less than the second element
cu_fun bool less(const uint32_t *element1, const size_t e1_size, const uint32_t *element2, const size_t e2_size)
{
    assert(e1_size == e2_size);
    for (size_t i = e2_size - 1; i > 0; i--)
        if (element1[i] > element2[i])
            return false;
        else if (element1[i] < element2[i])
            return true;
    return element1[0] < element2[0];
}

// Returns the carry, true if there was a carry, false otherwise
cu_fun bool _add(uint32_t *element1, const size_t e1_size, const uint32_t *element2, const size_t e2_size)
{
    assert(e1_size == e2_size);
    uint32_t carry = 0;
    for (int i = 0; i < e1_size; i++)
    {
        uint64_t tmp = (uint64_t)element1[i];
        tmp += carry;
        tmp += (uint64_t)element2[i];
        element1[i] = (uint32_t)(tmp);
        carry = (uint32_t)(tmp >> 32);
    }
    return carry;
}

// Fails if the second number is bigger than the first
cu_fun bool _subtract(uint32_t *element1, const size_t e1_size, bool carry, const uint32_t *element2, const size_t e2_size)
{
    assert(e1_size == e2_size);
    bool borrow = false;
    for (int i = 0; i < e1_size; i++)
    {
        uint64_t tmp = (uint64_t)element1[i];
        bool underflow = (tmp == 0) && (element2[i] > 0 || borrow);
        if (borrow)
            tmp--;
        borrow = underflow || (tmp < element2[i]);
        if (borrow)
            tmp += ((uint64_t)1 << 33);
        element1[i] = tmp - element2[i];
    }
    //assert(borrow == carry);
    return borrow;
}

cu_fun void montyNormalize(uint32_t *result, const size_t a_size, const uint32_t *mod, const bool msb)
{
    uint32_t u[SIZE] = {0};
    memcpy(u, result, a_size);
    bool borrow = _subtract(u, SIZE, false, mod, SIZE);
    if (msb || !borrow)
    {
        assert(!msb || msb == borrow);
        memcpy(result, u, a_size);
    }
}

cu_fun void ciosMontgomeryMultiply(uint32_t *result,
                                   const uint32_t *a, const size_t a_size,
                                   const uint32_t *b, const uint32_t *mod)
{
    uint64_t temp;
    for (size_t i = 0; i < a_size; i++)
    {
        uint32_t carry = 0;
        for (size_t j = 0; j < a_size; j++)
        {
            temp = result[j];
            temp += (uint64_t)a[j] * (uint64_t)b[i];
            temp += carry;
            result[j] = (uint32_t)temp;
            carry = temp >> 32;
        }
        temp = result[a_size] + carry;
        result[a_size] = (uint32_t)temp;
        result[a_size + 1] = temp >> 32;
        uint32_t m = (uint32_t)((uint64_t)result[0] * m_inv);
        temp = result[0] + (uint64_t)m * (uint64_t)mod[0];
        carry = temp >> 32;
        for (size_t j = 1; j < a_size; j++)
        {
            temp = result[j];
            temp += (uint64_t)m * (uint64_t)mod[j];
            temp += carry;
            result[j - 1] = (uint32_t)temp;
            carry = temp >> 32;
        }
        temp = result[a_size] + carry;
        result[a_size - 1] = (uint32_t)temp;
        result[a_size] = temp >> 32;
    }
    bool msb = result[a_size] > 0;
    montyNormalize(result, a_size, mod, msb);
}

//Adds two elements
cu_fun void Scalar::add(Scalar &fld1, const Scalar &fld2) const
{
    bool carry = _add(fld1.im_rep, SIZE, fld2.im_rep, SIZE);
    if (carry || less(_mod, SIZE, fld1.im_rep, SIZE))
        _subtract(fld1.im_rep, SIZE, false, _mod, SIZE);
}

//Subtract element two from element one
cu_fun void Scalar::subtract(Scalar &fld1, const Scalar &fld2) const
{
    bool carry = false;
    if (less(fld1.im_rep, SIZE, fld2.im_rep, SIZE))
        carry = _add(fld1.im_rep, SIZE, _mod, SIZE);
    _subtract(fld1.im_rep, SIZE, carry, fld2.im_rep, SIZE);
}


//Multiply two elements
cu_fun void Scalar::mul(Scalar &fld1, const Scalar &fld2, const uint32_t *mod) const
{

    uint32_t tmp[SIZE + 2] = {0};
    ciosMontgomeryMultiply(tmp, fld1.im_rep, SIZE, fld2.im_rep, mod);
    for (size_t i = 0; i < SIZE; i++)
        fld1.im_rep[i] = tmp[i];
}

cu_fun void to_monty(Scalar &a)
{
    // a = a << 2^(32*SIZE)
    // a = a % _mod
}

cu_fun void from_monty(Scalar &a)
{
    Scalar s = Scalar::one();
    a = a * s;
}

//Exponentiates this element
cu_fun void Scalar::pow(Scalar &fld1, const uint32_t pow) const
{
    if (pow == 0)
    {
        fld1 = Scalar::one();
        return;
    }

    if (pow == 1)
    {
        return;
    }

    uint32_t tmp[SIZE + 2];
    uint32_t temp[SIZE];

    to_monty(fld1);

    for (size_t i = 0; i < SIZE; i++)
        temp[i] = fld1.im_rep[i];

    for (size_t i = 0; i < pow - 1; i++)
    {
        memset(tmp, 0, (SIZE + 2) * sizeof(uint32_t));
        ciosMontgomeryMultiply(tmp, fld1.im_rep, SIZE, temp, _mod);
        for (size_t k = 0; k < SIZE; k++)
            fld1.im_rep[k] = tmp[k];
    }
    from_monty(fld1);
}
