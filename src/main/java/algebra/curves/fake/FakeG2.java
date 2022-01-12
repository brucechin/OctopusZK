/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.fake;

import algebra.fields.Fp;
import algebra.curves.AbstractG2;
import algebra.curves.fake.abstract_fake_parameters.AbstractFakeG2Parameters;

import java.io.Serializable;
import java.util.ArrayList;
import algebra.fields.Fp2;
import java.math.BigInteger;

public class FakeG2 extends AbstractG2<FakeG2> implements Serializable {

    public Fp element;
    public AbstractFakeG2Parameters FakeG2Parameters;

    public FakeG2(final Fp element, final AbstractFakeG2Parameters FakeG2Parameters) {
        this.element = element;
        this.FakeG2Parameters = FakeG2Parameters;
    }

    public FakeG2(final BigInteger number, final AbstractFakeG2Parameters FakeG2Parameters) {
        this(new Fp(number, FakeG2Parameters.FqParameters()), FakeG2Parameters);
    }

    public FakeG2(final String number, final AbstractFakeG2Parameters FakeG2Parameters) {
        this(new Fp(number, FakeG2Parameters.FqParameters()), FakeG2Parameters);
    }

    public FakeG2(final long number, final AbstractFakeG2Parameters FakeG2Parameters) {
        this(new Fp(Long.toString(number), FakeG2Parameters.FqParameters()), FakeG2Parameters);
    }

    public FakeG2 self() {
        return this;
    }

    public FakeG2 add(final FakeG2 that) {
        return new FakeG2(this.element.add(that.element), FakeG2Parameters);
    }

    public FakeG2 sub(final FakeG2 that) {
        return new FakeG2(this.element.sub(that.element), FakeG2Parameters);
    }

    public FakeG2 zero() {
        return FakeG2Parameters.ZERO();
    }

    public boolean isZero() {
        return this.equals(FakeG2Parameters.ZERO());
    }

    public boolean isSpecial() {
        return isZero();
    }

    public FakeG2 one() {
        return FakeG2Parameters.ONE();
    }

    public boolean isOne() {
        return this.equals(FakeG2Parameters.ONE());
    }

    public FakeG2 random(final Long seed, final byte[] secureSeed) {
        return new FakeG2(element.random(seed, secureSeed), FakeG2Parameters);
    }

    public FakeG2 negate() {
        return new FakeG2(this.element.negate(), FakeG2Parameters);
    }

    public FakeG2 twice() {
        return this.add(this);
    }

    public int bitSize() {
        return this.element.bitSize();
    }

    public ArrayList<Integer> fixedBaseWindowTable() {
        return FakeG2Parameters.fixedBaseWindowTable();
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
    public  ArrayList<BigInteger> BN254G2ToBigInteger(){
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        return result;
    }

    public boolean equals(final FakeG2 that) {
        if (that == null) {
            return false;
        }

        return this.element.equals(that.element);
    }

}
