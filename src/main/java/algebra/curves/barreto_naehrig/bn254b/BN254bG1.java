/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.barreto_naehrig.bn254b;

import algebra.curves.barreto_naehrig.BNG1;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields.BN254bFq;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields.BN254bFr;
import algebra.curves.barreto_naehrig.bn254b.bn254b_parameters.BN254bG1Parameters;

import java.util.ArrayList;
import algebra.fields.Fp2;
import algebra.math.BigInteger;
public class BN254bG1
        extends BNG1<BN254bFr, BN254bFq, BN254bG1, BN254bG1Parameters> {

    public static final BN254bG1Parameters G1Parameters = new BN254bG1Parameters();

    public BN254bG1(final BN254bFq X, final BN254bFq Y, final BN254bFq Z) {
        super(X, Y, Z, G1Parameters);
    }

    public BN254bG1 self() {
        return this;
    }

    public BN254bG1 construct(final BN254bFq X, final BN254bFq Y, final BN254bFq Z) {
        return new BN254bG1(X, Y, Z);
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

