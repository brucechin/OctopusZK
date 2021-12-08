/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.fake;

import algebra.fields.Fp;
import algebra.curves.AbstractPairing;

import java.math.BigInteger;

import static algebra.curves.fake.FakeInitialize.GTParameters;

public class FakePairing extends AbstractPairing<FakeG1, FakeG2, FakeGT> {

    public FakeG1Precompute precomputeG1(final FakeG1 P) {
        return new FakeG1Precompute(P);
    }

    public FakeG2Precompute precomputeG2(final FakeG2 Q) {
        return new FakeG2Precompute(Q);
    }

    private Fp millerLoop(final FakeG1Precompute PPrec, final FakeG2Precompute QPrec) {
        final BigInteger PQ = PPrec.P.toBigInteger().multiply(QPrec.Q.toBigInteger());
        return new Fp(PQ, GTParameters.FqParameters());
    }

    private Fp atePairing(final FakeG1 P, final FakeG2 Q) {
        final FakeG1Precompute PPrec = precomputeG1(P);
        final FakeG2Precompute QPrec = precomputeG2(Q);
        return millerLoop(PPrec, QPrec);
    }

    public FakeGT reducedPairing(final FakeG1 P, final FakeG2 Q) {
        final Fp f = atePairing(P, Q);
        return new FakeGT(f.toBigInteger(), GTParameters);
    }

    public class FakeG1Precompute extends G1Precompute {
        FakeG1 P;

        private FakeG1Precompute(final FakeG1 _P) {
            P = _P;
        }

        public boolean equals(final FakeG1Precompute that) {
            return P.equals(that.P);
        }
    }

    public class FakeG2Precompute extends G2Precompute {
        FakeG2 Q;

        private FakeG2Precompute(final FakeG2 _Q) {
            Q = _Q;
        }

        public boolean equals(final FakeG2Precompute that) {
            return Q.equals(that.Q);
        }
    }
}
