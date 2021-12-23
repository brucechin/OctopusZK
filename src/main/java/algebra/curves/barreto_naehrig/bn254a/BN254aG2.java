/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.barreto_naehrig.bn254a;

import algebra.curves.barreto_naehrig.BNG2;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFq;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFq2;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFr;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG2Parameters;
import java.math.BigInteger;
import java.util.ArrayList;

public class BN254aG2 extends BNG2<BN254aFr, BN254aFq, BN254aFq2, BN254aG2, BN254aG2Parameters> {

    private static final BN254aG2Parameters G2Parameters = new BN254aG2Parameters();

    public BN254aG2(
            final BN254aFq2 X,
            final BN254aFq2 Y,
            final BN254aFq2 Z) {
        super(X, Y, Z, G2Parameters);
    }

    public BN254aG2 self() {
        return this;
    }

    public BN254aG2 construct(final BN254aFq2 X, final BN254aFq2 Y, final BN254aFq2 Z) {
        return new BN254aG2(X, Y, Z);
    }

    public BigInteger toBigInteger() {
        //WATNING do not call this function
        return BigInteger.ZERO;
    }

    public void setBigInteger(BigInteger bigInteger) {
        //WATNING do not call this function
    }


    public ArrayList<BigInteger> BN254G1ToBigInteger() {
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        //WARNING: do not use this for BN254G2
        return result;
    }


    public void setBigIntegerBN254G1(BigInteger x, BigInteger y, BigInteger z) {
        //WARNING: do not use this for BN254G2
    }

    public void setBigIntegerBN254G2(BigInteger x1, BigInteger x2, BigInteger y1, BigInteger y2, BigInteger z1, BigInteger z2) {
        this.X.element().c0.number = x1;
        this.X.element().c1.number = x2;
        this.Y.element().c0.number = y1;
        this.Y.element().c1.number = y2;        
        this.Z.element().c0.number = z1;
        this.Z.element().c1.number = z2;
    }

    public ArrayList<BigInteger> BN254G2ToBigInteger() {
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        result.add(this.X.element().c0.toBigInteger());
        result.add(this.X.element().c1.toBigInteger());
        result.add(this.Y.element().c0.toBigInteger());
        result.add(this.Y.element().c1.toBigInteger());
        result.add(this.Z.element().c0.toBigInteger());
        result.add(this.Z.element().c1.toBigInteger());
        return result;
    }
}
