
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <vector>
#include <string>
#include <bitset>
#include "BigInt.h"
#include "Fp.h"


class AbstractFpParameters {
public:
    virtual BigInt modulus()= 0;

    virtual BigInt root()= 0;

    virtual Fp multiplicativeGenerator()= 0;

    virtual long numBits()= 0;

    virtual BigInt euler()= 0;

    virtual long s()= 0;

    virtual BigInt t()= 0;

    virtual BigInt tMinus1Over2()= 0;

    virtual Fp nqr()= 0;

    virtual Fp nqrTot()= 0;

    virtual Fp ZERO()= 0;

    virtual Fp ONE()= 0;
}


class BN254aFqParameters : AbstractFpParameters{
    BigInt modulus;
    BigInt root;
    Fp multiplicativeGenerator;
    long numBits;

    BigInt euler;
    long s;
    BigInt t;
    BigInt tMinus1Over2;
    Fp nqr;
    Fp nqrTot;

    Fp ZERO;
    Fp ONE;

    BN254aFqParameters() {
        this.modulus =  BigInt("21888242871839275222246405745257275088696311157297823662689037894645226208583");
        this.root =  BigInt("21888242871839275222246405745257275088696311157297823662689037894645226208582");
        this.multiplicativeGenerator =  Fp("3", this);
        this.numBits = 254;

        this.euler =  BigInt("10944121435919637611123202872628637544348155578648911831344518947322613104291");
        this.s = 1;
        this.t =  BigInt("10944121435919637611123202872628637544348155578648911831344518947322613104291");
        this.tMinus1Over2 =  BigInt("5472060717959818805561601436314318772174077789324455915672259473661306552145");
        this.nqr =  Fp("3", this);
        this.nqrTot =  Fp("21888242871839275222246405745257275088696311157297823662689037894645226208582", this);

        this.ZERO =  Fp(BigInt.ZERO(), this);
        this.ONE =  Fp(BigInt.ONE(), this);
    }

     BigInt modulus() {
        return modulus;
    }

     BigInt root() {
        return root;
    }

     Fp multiplicativeGenerator() {
        return multiplicativeGenerator;
    }

     long numBits() {
        return numBits;
    }

     BigInt euler() {
        return euler;
    }

     long s() {
        return s;
    }

     BigInt t() {
        return t;
    }

     BigInt tMinus1Over2() {
        return tMinus1Over2;
    }

     Fp nqr() {
        return nqr;
    }

     Fp nqrTot() {
        return nqrTot;
    }

     Fp ZERO() {
        return ZERO;
    }

     Fp ONE() {
        return ONE;
    }
};




class BN254aFrParameters extends AbstractBNFrParameters implements Serializable {
    BigInt modulus;
    BigInt root;
    Fp multiplicativeGenerator;
    long numBits;

    BigInt euler;
    long s;
    BigInt t;
    BigInt tMinus1Over2;
    Fp nqr;
    Fp nqrTot;

    Fp ZERO;
    Fp ONE;

    BN254aFrParameters() {
        this.modulus = new BigInt("21888242871839275222246405745257275088548364400416034343698204186575808495617");
        this.root = new BigInt("19103219067921713944291392827692070036145651957329286315305642004821462161904");
        this.multiplicativeGenerator = new Fp("5", this);
        this.numBits = 254;

        this.euler = new BigInt("10944121435919637611123202872628637544274182200208017171849102093287904247808");
        this.s = 28;
        this.t = new BigInt("81540058820840996586704275553141814055101440848469862132140264610111");
        this.tMinus1Over2 = new BigInt("40770029410420498293352137776570907027550720424234931066070132305055");
        this.nqr = new Fp("5", this);
        this.nqrTot = new Fp("19103219067921713944291392827692070036145651957329286315305642004821462161904", this);

        this.ZERO = new Fp(BigInt.ZERO(), this);
        this.ONE = new Fp(BigInt.ONE(), this);
    }

    BigInt modulus() {
        return modulus;
    }

    BigInt root() {
        return root;
    }

    Fp multiplicativeGenerator() {
        return multiplicativeGenerator;
    }

    long numBits() {
        return numBits;
    }

    BigInt euler() {
        return euler;
    }

    long s() {
        return s;
    }

    BigInt t() {
        return t;
    }

    BigInt tMinus1Over2() {
        return tMinus1Over2;
    }

    Fp nqr() {
        return nqr;
    }

    Fp nqrTot() {
        return nqrTot;
    }

    Fp ZERO() {
        return ZERO;
    }

    Fp ONE() {
        return ONE;
    }
}


 /* Scalar field Fr */
 class BN254aFr{

        static const BN254aFrParameters FrParameters = new BN254aFrParameters();
        static const BN254aFr ZERO = new BN254aFr(FrParameters.ZERO());
        static const BN254aFr ONE = new BN254aFr(FrParameters.ONE());
        static const BN254aFr MULTIPLICATIVE_GENERATOR =
                new BN254aFr(FrParameters.multiplicativeGenerator());

        Fp element;

        BN254aFr(const BigInt number) {
            this.element = new Fp(number, FrParameters);
        }

        BN254aFr(const Fp number) {
            this(number.toBigInteger());
        }

        BN254aFr(const String number) {
            this(new BigInt(number));
        }

        BN254aFr(const long number) {
            this(BigInt.valueOf(number));
        }

        BN254aFr self() {
            return this;
        }

        Fp element() {
            return element;
        }

        BN254aFr zero() {
            return ZERO;
        }

        BN254aFr one() {
            return ONE;
        }

        BN254aFr multiplicativeGenerator() {
            return MULTIPLICATIVE_GENERATOR;
        }

        BN254aFr construct(const long number) {
            return new BN254aFr(number);
        }

        BN254aFr construct(const Fp element) {
            return new BN254aFr(element);
        }

    }

 /* Base field Fq */
class BN254aFq {

        static const BN254aFqParameters FqParameters = new BN254aFqParameters();
        static const BN254aFq ZERO = new BN254aFq(FqParameters.ZERO());
        static const BN254aFq ONE = new BN254aFq(FqParameters.ONE());
        static const BN254aFq MULTIPLICATIVE_GENERATOR =
                new BN254aFq(FqParameters.multiplicativeGenerator());

        Fp element;

        BN254aFq(const Fp element) {
            element = element;
        }

        BN254aFq(const BigInt number) {
            element = new Fp(number, FqParameters);
        }

        BN254aFq(const string number) {
            this(new BigInt(number));
        }

        Fp element() {
            return element;
        }

        BN254aFq zero() {
            return ZERO;
        }

        BN254aFq one() {
            return ONE;
        }

        BN254aFq multiplicativeGenerator() {
            return MULTIPLICATIVE_GENERATOR;
        }

        BN254aFq construct(const Fp element) {
            return new BN254aFq(element);
        }

         BN254aFq construct(const string element) {
            return new BN254aFq(element);
        }

        BN254aFq construct(const long number) {
            return new BN254aFq(number);
        }

         String toString() {
            return this.element.toString();
        }
    };