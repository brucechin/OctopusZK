
#include<iostream>
#include<stdexcept>
#include<unistd.h>
#include<cstring>
#include <vector>
#include <string>
#include <bitset>
#include "BigInt.h"
#include "Fp.h"
#include "BN254aFields.h"


 class BNG1{
    const BNG1ParametersT G1Parameters;
    const BNFqT X;
    const BNFqT Y;
    const BNFqT Z;

    BNG1(const BNFqT X, const BNFqT Y, const BNFqT Z, const BNG1ParametersT G1Parameters) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.G1Parameters = G1Parameters;
    }

    abstract BNG1T construct(BNFqT X, BNFqT Y, BNFqT Z);

    BNG1T add(const BNG1T that) {
        assert (that != null);

        // Handle special cases having to do with O
        if (isZero()) {
            return that;
        }

        if (that.isZero()) {
            return this.self();
        }

        // No need to handle points of order 2,4
        // (they cannot exist in a modulus-order subgroup)

        // Check for doubling case

        // Using Jacobian coordinates so:
        // (X1:Y1:Z1) = (X2:Y2:Z2)
        // iff
        // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
        // iff
        // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

        const BNFqT Z1Z1 = this.Z.square();
        const BNFqT Z2Z2 = that.Z.square();

        const BNFqT U1 = this.X.mul(Z2Z2);
        const BNFqT U2 = that.X.mul(Z1Z1);

        const BNFqT Z1_cubed = this.Z.mul(Z1Z1);
        const BNFqT Z2_cubed = that.Z.mul(Z2Z2);

        const BNFqT S1 = this.Y.mul(Z2_cubed);      // S1 = Y1 * Z2 * Z2Z2
        const BNFqT S2 = that.Y.mul(Z1_cubed);      // S2 = Y2 * Z1 * Z1Z1

        if (U1.equals(U2) && S1.equals(S2)) {
            // Double case; nothing above can be reused.
            return twice();
        }

        // Rest of the add case.
        const BNFqT H = U2.sub(U1);                                   // H = U2-U1
        const BNFqT S2_minus_S1 = S2.sub(S1);
        const BNFqT I = H.add(H).square();                            // I = (2 * H)^2
        const BNFqT J = H.mul(I);                                     // J = H * I
        const BNFqT r = S2_minus_S1.add(S2_minus_S1);                 // r = 2 * (S2-S1)
        const BNFqT V = U1.mul(I);                                    // V = U1 * I
        const BNFqT X3 = r.square().sub(J).sub(V.add(V));             // X3 = r^2 - J - 2 * V
        const BNFqT S1_J = S1.mul(J);
        const BNFqT Y3 = r.mul(V.sub(X3)).sub(S1_J.add(S1_J));        // Y3 = r * (V-X3)-2 * S1_J
        const BNFqT Z3 = this.Z.add(that.Z).square().sub(Z1Z1).sub(Z2Z2)
                .mul(H);                                                  // Z3 = ((Z1+Z2)^2-Z1Z1-Z2Z2) * H

        return this.construct(X3, Y3, Z3);
    }

    BNG1T sub(const BNG1T that) {
        return this.add(that.negate());
    }

    bool isZero() {
        return this.Z.isZero();
    }

    bool isSpecial() {
        return isZero() || isOne();
    }

    bool isOne() {
        return this.X.equals(this.one().X)
                && this.Y.equals(this.one().Y)
                && this.Z.equals(this.one().Z);
    }

    BNG1T zero() {
        return this.G1Parameters.ZERO();
    }

    BNG1T one() {
        return this.G1Parameters.ONE();
    }

    BNG1T random(const Long seed, const byte[] secureSeed) {
        return this.one().mul(this.G1Parameters.oneFr().random(seed, secureSeed));
    }

    BNG1T negate() {
        return this.construct(this.X, this.Y.negate(), this.Z);
    }

    BNG1T twice() {
        // handle point at infinity
        if (isZero()) {
            return this.self();
        }

        // No need to handle points of order 2,4
        // (they cannot exist in a modulus-order subgroup)

        // NOTE: does not handle O and pts of order 2,4
        // http://www.hyperelliptic.org/EFD/g1p/auto-shortw-jacobian-0.html#doubling-dbl-2009-l

        const BNFqT A = this.X.square();                       // A = X1^2
        const BNFqT B = this.Y.square();                       // B = Y1^2
        const BNFqT C = B.square();                            // C = B^2
        BNFqT D = this.X.add(B).square().sub(A).sub(C);
        D = D.add(D);                                          // D = 2 * ((X1 + B)^2 - A - C)
        const BNFqT E = A.add(A).add(A);                       // E = 3 * A
        const BNFqT F = E.square();                            // F = E^2
        const BNFqT X3 = F.sub(D.add(D));                      // X3 = F - 2 D
        BNFqT eightC = C.add(C);
        eightC = eightC.add(eightC);
        eightC = eightC.add(eightC);
        const BNFqT Y3 = E.mul(D.sub(X3)).sub(eightC);         // Y3 = E * (D - X3) - 8 * C
        const BNFqT Y1Z1 = this.Y.mul(this.Z);
        const BNFqT Z3 = Y1Z1.add(Y1Z1);                       // Z3 = 2 * Y1 * Z1

        return this.construct(X3, Y3, Z3);
    }

    BNG1T toAffineCoordinates() {
        if (isZero()) {
            return this.construct(this.X.zero(), this.Y.one(), this.Z.zero());
        } else {
            BNFqT ZInverse = this.Z.inverse();
            BNFqT Z2Inverse = ZInverse.square();
            BNFqT Z3Inverse = Z2Inverse.mul(ZInverse);
            return this.construct(this.X.mul(Z2Inverse), this.Y.mul(Z3Inverse), this.Z.one());
        }
    }

    int bitSize() {
        return Math.max(this.X.bitSize(), Math.max(this.Y.bitSize(), this.Z.bitSize()));
    }

    ArrayList<Integer> fixedBaseWindowTable() {
        return this.G1Parameters.fixedBaseWindowTable();
    }


    bool equals(const BNG1T that) {
        if (isZero()) {
            return that.isZero();
        }

        if (that.isZero()) {
            return false;
        }

        // Now neither is O.

        // using Jacobian coordinates so:
        // (X1:Y1:Z1) = (X2:Y2:Z2)
        // iff
        // X1/Z1^2 == X2/Z2^2 and Y1/Z1^3 == Y2/Z2^3
        // iff
        // X1 * Z2^2 == X2 * Z1^2 and Y1 * Z2^3 == Y2 * Z1^3

        const BNFqT Z1_squared = this.Z.square();
        const BNFqT Z2_squared = that.Z.square();

        if (!this.X.mul(Z2_squared).equals(that.X.mul(Z1_squared))) {
            return false;
        }

        const BNFqT Z1_cubed = this.Z.mul(Z1_squared);
        const BNFqT Z2_cubed = that.Z.mul(Z2_squared);

        if (!this.Y.mul(Z2_cubed).equals(that.Y.mul(Z1_cubed))) {
            return false;
        }

        return true;
    }
}
