/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.groups;

import algebra.fields.AbstractFieldElementExpanded;

import java.io.Serializable;
import java.util.ArrayList;
import algebra.fields.Fp2;
import java.math.BigInteger;

public abstract class AbstractGroup<GroupT extends AbstractGroup<GroupT>> implements Serializable {

    /* Returns self element */
    public abstract GroupT self();

    /* Returns this + that */
    public abstract GroupT add(final GroupT that);

    /* Returns this - that */
    public abstract GroupT sub(final GroupT that);

    /* Returns (BigInteger) scalar * this */
    public GroupT mul(final BigInteger scalar) {
        GroupT base = this.self();

        if (scalar.equals(BigInteger.ONE)) {
            return base;
        }

        GroupT result = base.zero();

        boolean found = false;
        for (int i = scalar.bitLength() - 1; i >= 0; i--) {
            if (found) {
                result = result.twice();
            }

            if (scalar.testBit(i)) {
                found = true;
                result = result.add(base);
            }
        }

        return result;
    }

    /* Returns (FieldT) scalar * this */
    public GroupT mul(final AbstractFieldElementExpanded that) {
        return this.mul(that.toBigInteger());
    }

    /* Returns the zero element of the group */
    public abstract GroupT zero();

    /* Returns the one element of the group */
    public abstract GroupT one();

    /* Returns this == zero */
    public abstract boolean isZero();

    /* Returns this == one */
    public abstract boolean isOne();

    /* Returns -this */
    public abstract GroupT negate();

    /* Returns this+this */
    public abstract GroupT twice();

    /* Fixed base window table for G1 and G2. */
    public abstract ArrayList<Integer> fixedBaseWindowTable();

    public abstract BigInteger toBigInteger();
    public abstract ArrayList<BigInteger> BN254G1ToBigInteger();
    public abstract ArrayList<BigInteger> BN254G2ToBigInteger();


    public abstract void setBigInteger(final BigInteger bigInteger);
    public abstract void setBigIntegerBN254G1(BigInteger x, BigInteger y, BigInteger z);
    public abstract void setBigIntegerBN254G2(BigInteger x1, BigInteger x2, BigInteger y1, BigInteger y2, BigInteger z1, BigInteger z2);
    public abstract ArrayList<Fp2> BN254G2ToFp2();
    /**
     * If secureSeed is provided, returns cryptographically secure random group element using byte[].
     * Else if seed is provided, returns pseudorandom group element using long as seed.
     * Else, returns a pseudorandom group element without a seed.
     */
    public abstract GroupT random(final Long seed, final byte[] secureSeed);

    /* Returns this == that */
    public abstract boolean equals(final GroupT that);

}
