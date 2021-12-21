/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.barreto_naehrig.bn254a;
import java.math.BigInteger;
import java.util.ArrayList;

import algebra.curves.barreto_naehrig.BNG1;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFq;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFr;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG1Parameters;

public class BN254aG1 extends BNG1<BN254aFr, BN254aFq, BN254aG1, BN254aG1Parameters> {

    public static final BN254aG1Parameters G1Parameters = new BN254aG1Parameters();

    public BN254aG1(
            final BN254aFq X,
            final BN254aFq Y,
            final BN254aFq Z) {
        super(X, Y, Z, G1Parameters);
    }

    public BN254aG1 self() {
        return this;
    }

    public BN254aG1 construct(final BN254aFq X, final BN254aFq Y, final BN254aFq Z) {
        return new BN254aG1(X, Y, Z);
    }

    public BigInteger toBigInteger() {
        return BigInteger.ZERO;
    }

    public ArrayList<BigInteger> BN254G1ToBigInteger() {
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        result.add(this.X.toBigInteger());
        result.add(this.Y.toBigInteger());
        result.add(this.Z.toBigInteger());
        return result;
    }

    public void setBigInteger(BigInteger bigInteger) {
        //do not use it
    }


    public void setBigIntegerBN254G1(BigInteger x, BigInteger y, BigInteger z) {
        this.X.element().number = x;
        this.Y.element().number = y;
        this.Z.element().number = z;
    }
}
