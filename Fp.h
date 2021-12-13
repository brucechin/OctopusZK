#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <vector>
#include <string>
#include <bitset>
#include "BigInt.h"
#include "BN254aFields.h"
class Fp {
    BigInt number;
    AbstractFpParameters FpParameters;

    Fp(const BigInt number, const AbstractFpParameters FpParameters) {
        this.number = number.mod(FpParameters.modulus());
        this.FpParameters = FpParameters;
    }

    Fp(const string number, const AbstractFpParameters FpParameters) {
        this(BigInt(number), FpParameters);
    }


    Fp add(const Fp that) {
        return new Fp(number +that.number, FpParameters);
    }

    Fp sub(const Fp that) {
        return new Fp(number.subtract(that.number), FpParameters);
    }

    Fp mul(const Fp that) {
        return new Fp(number * that.number), FpParameters);
    }

    Fp mul(const BigInt that) {
        return new Fp(number.multiply(that), FpParameters);
    }

    Fp zero() {
        return FpParameters.ZERO();
    }

    bool isZero() {
        return equals(FpParameters.ZERO());
    }

    Fp one() {
        return FpParameters.ONE();
    }

    bool isOne() {
        return equals(FpParameters.ONE());
    }

    // Fp negate() {
    //     return  Fp(number.negate(), FpParameters);
    // }

    Fp square() {
        return  Fp(number * number, FpParameters);
    }

    Fp inverse() {
        return new Fp(number.modInverse(FpParameters.modulus()), FpParameters);
    }

    Fp multiplicativeGenerator() {
        return FpParameters.multiplicativeGenerator();
    }

    Fp rootOfUnity(const long order) {
        const BigInt exponent = FpParameters.modulus().divide(BigInt.valueOf(order));

        return  Fp(FpParameters.root().modPow(exponent, FpParameters.modulus()), FpParameters);
    }

    // int bitSize() {
    //     return number.bitLength();
    // }

    Fp construct(const long value) {
        return  Fp(value, FpParameters);
    }

    BigInt toBigInteger() {
        return number;
    }

    void setBigInteger(BigInt bigInteger) {
        this.number = bigInteger;
    }

    boolean equals(const Fp that) {
        if (that == null) {
            return false;
        }

        return number.equals(that.number);
    }
}
