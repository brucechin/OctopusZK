/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.barreto_naehrig.bn254b;

import algebra.curves.barreto_naehrig.BNG2;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields.BN254bFq;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields.BN254bFq2;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields.BN254bFr;
import algebra.curves.barreto_naehrig.bn254b.bn254b_parameters.BN254bG2Parameters;

import java.util.ArrayList;
import algebra.fields.Fp2;
import java.math.BigInteger;
public class BN254bG2
        extends BNG2<BN254bFr, BN254bFq, BN254bFq2, BN254bG2, BN254bG2Parameters> {

    private static final BN254bG2Parameters G2Parameters = new BN254bG2Parameters();

    public BN254bG2(final BN254bFq2 X, final BN254bFq2 Y, final BN254bFq2 Z) {
        super(X, Y, Z, G2Parameters);
    }

    public BN254bG2 self() {
        return this;
    }

    public BN254bG2 construct(final BN254bFq2 X, final BN254bFq2 Y, final BN254bFq2 Z) {
        return new BN254bG2(X, Y, Z);
    }

    public BigInteger toBigInteger() {
        return BigInteger.ZERO;
    }

    public void setBigInteger(BigInteger bigInteger) {
        //do not use it
    }

    public  ArrayList<Fp2> BN254G2ToFp2(){
        ArrayList<Fp2> res = new ArrayList<>();
        return res;
    }

    public  void setBigIntegerBN254G2(BigInteger x1, BigInteger x2, BigInteger y1, BigInteger y2, BigInteger z1, BigInteger z2){

    }


    public  void setBigIntegerBN254G1(BigInteger x, BigInteger y, BigInteger z){

    }

    public  ArrayList<BigInteger> BN254G1ToBigInteger(){
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        return result;
    }
    public  ArrayList<BigInteger> BN254G2ToBigInteger(){
        ArrayList<BigInteger> result = new ArrayList<BigInteger>();
        return result;
    }

}
