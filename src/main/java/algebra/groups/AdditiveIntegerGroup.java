/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.groups;

import algebra.groups.abstractintegergroupparameters.AbstractAdditiveIntegerGroupParameters;

import java.math.BigInteger;
import java.security.SecureRandom;
import java.util.ArrayList;
import java.util.Random;

import algebra.fields.Fp2;

/*
 * AdditiveIntegerGroup is an additive cyclic group of order r.
 * If the order is not provided for the constructor, the group behaves as the integer group Z.
 */
public class AdditiveIntegerGroup extends AbstractGroup<AdditiveIntegerGroup> {
    protected BigInteger number;
    private final AbstractAdditiveIntegerGroupParameters GroupParameters;
    private final ArrayList<Integer> fixedBaseWindowTable = new ArrayList<>(0);

    public AdditiveIntegerGroup(
            final BigInteger number,
            final AbstractAdditiveIntegerGroupParameters GroupParameters) {
        if (GroupParameters.modulus() == null) {
            this.number = number;
        } else {
            this.number = number.mod(GroupParameters.modulus());
        }

        this.GroupParameters = GroupParameters;
    }

    public AdditiveIntegerGroup(
            final String number,
            final AbstractAdditiveIntegerGroupParameters GroupParameters) {
        this(new BigInteger(number), GroupParameters);
    }

    public AdditiveIntegerGroup(
            final long number,
            final AbstractAdditiveIntegerGroupParameters GroupParameters) {
        this(new BigInteger(Long.toString(number)), GroupParameters);
    }

    public AdditiveIntegerGroup self() {
        return this;
    }

    public AdditiveIntegerGroup add(final AdditiveIntegerGroup that) {
        return new AdditiveIntegerGroup(this.number.add(that.number), GroupParameters);
    }

    public AdditiveIntegerGroup sub(final AdditiveIntegerGroup that) {
        return new AdditiveIntegerGroup(this.number.subtract(that.number), GroupParameters);
    }

    public AdditiveIntegerGroup mul(final BigInteger scalar) {
        if (scalar.equals(BigInteger.ONE)) {
            return this;
        }

        AdditiveIntegerGroup result = zero();

        boolean found = false;
        for (int i = scalar.bitLength() - 1; i >= 0; i--) {
            if (found) {
                result = result.twice();
            }

            if (scalar.testBit(i)) {
                found = true;
                result = result.add(this);
            }
        }

        return result;
    }

    public void setBigInteger(BigInteger bigInteger) {
        this.number = bigInteger;
    }


    public AdditiveIntegerGroup zero() {
        return GroupParameters.ZERO();
    }

    public AdditiveIntegerGroup one() {
        return GroupParameters.ONE();
    }

    public boolean isZero() {
        return this.equals(GroupParameters.ZERO());
    }

    public boolean isOne() {
        return this.equals(GroupParameters.ONE());
    }

    public AdditiveIntegerGroup random(final Long seed, final byte[] secureSeed) {
        if (secureSeed != null && secureSeed.length > 0) {
            return new AdditiveIntegerGroup(new SecureRandom(secureSeed).nextLong(), GroupParameters);
        } else if (seed != null) {
            return new AdditiveIntegerGroup(new Random(seed).nextLong(), GroupParameters);
        } else {
            return new AdditiveIntegerGroup(new Random().nextLong(), GroupParameters);
        }
    }

    public AdditiveIntegerGroup negate() {
        if (GroupParameters.modulus() == null) {
            return new AdditiveIntegerGroup(this.number.negate(), GroupParameters);
        }
        return new AdditiveIntegerGroup(GroupParameters.modulus().subtract(number), GroupParameters);
    }

    public AdditiveIntegerGroup twice() {
        return new AdditiveIntegerGroup(this.number.add(number), GroupParameters);
    }

    public ArrayList<Integer> fixedBaseWindowTable() {
        return fixedBaseWindowTable;
    }

    public String toString() {
        return this.number.toString();
    }

    public BigInteger toBigInteger() {
        return this.number;
    }

    public boolean equals(final AdditiveIntegerGroup that) {
        if (that == null) {
            return false;
        }
        return this.number.equals(that.number);
    }


    public  ArrayList<Fp2> BN254G2ToFp2(){
        ArrayList<Fp2> res = new ArrayList<>();
        return res;
    }
    public  void setBigIntegerBN254G1(BigInteger x, BigInteger y, BigInteger z){

    }
    public  void setBigIntegerBN254G2(BigInteger x1, BigInteger x2, BigInteger y1, BigInteger y2, BigInteger z1, BigInteger z2){

    }

    public  ArrayList<BigInteger> BN254G1ToBigInteger(){
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        return result;
    }
    public ArrayList<BigInteger> BN254G2ToBigInteger(){
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        return result;
    }
    
}
