/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.fake;

import algebra.fields.Fp;
import algebra.curves.AbstractG1;
import algebra.curves.fake.abstract_fake_parameters.AbstractFakeG1Parameters;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;
import algebra.fields.Fp2;

public class FakeG1 extends AbstractG1<FakeG1> implements Serializable {

    public  Fp element;
    public AbstractFakeG1Parameters FakeG1Parameters;

    public FakeG1(final Fp element, final AbstractFakeG1Parameters FakeG1Parameters) {
        this.element = element;
        this.FakeG1Parameters = FakeG1Parameters;
    }

    public FakeG1(final BigInteger number, final AbstractFakeG1Parameters FakeG1Parameters) {
        this(new Fp(number, FakeG1Parameters.FqParameters()), FakeG1Parameters);
    }

    public FakeG1(final String number, final AbstractFakeG1Parameters FakeG1Parameters) {
        this(new Fp(number, FakeG1Parameters.FqParameters()), FakeG1Parameters);
    }

    public FakeG1(final long number, final AbstractFakeG1Parameters FakeG1Parameters) {
        this(new Fp(Long.toString(number), FakeG1Parameters.FqParameters()), FakeG1Parameters);
    }

    public FakeG1 self() {
        return this;
    }

    public FakeG1 add(final FakeG1 that) {
        return new FakeG1(this.element.add(that.element), FakeG1Parameters);
    }

    public FakeG1 sub(final FakeG1 that) {
        return new FakeG1(this.element.sub(that.element), FakeG1Parameters);
    }

    public FakeG1 zero() {
        return FakeG1Parameters.ZERO();
    }

    public boolean isZero() {
        return this.equals(FakeG1Parameters.ZERO());
    }

    public boolean isSpecial() {
        return isZero();
    }

    public FakeG1 one() {
        return FakeG1Parameters.ONE();
    }

    public boolean isOne() {
        return this.equals(FakeG1Parameters.ONE());
    }

    public FakeG1 random(final Long seed, final byte[] secureSeed) {
        return new FakeG1(element.random(seed, secureSeed), FakeG1Parameters);
    }

    public FakeG1 negate() {
        return new FakeG1(this.element.negate(), FakeG1Parameters);
    }

    public FakeG1 twice() {
        return this.add(this);
    }

    public int bitSize() {
        return this.element.bitSize();
    }

    public ArrayList<Integer> fixedBaseWindowTable() {
        return FakeG1Parameters.fixedBaseWindowTable();
    }

    public BigInteger toBigInteger() {
        return this.element.toBigInteger();
    }

    public void setBigInteger(BigInteger bigInteger) {
        this.element.number = bigInteger;
    }

    public String toString() {
        return this.element.toString();
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

    public boolean equals(final FakeG1 that) {
        if (that == null) {
            return false;
        }

        return this.element.equals(that.element);
    }

}
