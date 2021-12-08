/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.fields;

import algebra.fields.abstractfieldparameters.AbstractFp6_2Over3_Parameters;

import java.io.Serializable;

public class Fp6_2Over3 extends AbstractFieldElement<Fp6_2Over3> implements Serializable {
    protected final Fp3 c0;
    protected final Fp3 c1;
    private final AbstractFp6_2Over3_Parameters Fp6Parameters;

    public Fp6_2Over3(final Fp3 c0, final Fp3 c1, final AbstractFp6_2Over3_Parameters Fp6Parameters) {
        this.c0 = c0;
        this.c1 = c1;
        this.Fp6Parameters = Fp6Parameters;
    }

    public Fp6_2Over3 self() {
        return this;
    }

    public Fp6_2Over3 add(final Fp6_2Over3 that) {
        return new Fp6_2Over3(c0.add(that.c0), c1.add(that.c1), Fp6Parameters);
    }

    public Fp6_2Over3 sub(final Fp6_2Over3 that) {
        return new Fp6_2Over3(c0.sub(that.c0), c1.sub(that.c1), Fp6Parameters);
    }

    public Fp6_2Over3 mul(final Fp that) {
        return new Fp6_2Over3(c0.mul(that), c1.mul(that), Fp6Parameters);
    }

    public Fp6_2Over3 mul(final Fp3 that) {
        return new Fp6_2Over3(c0.mul(that), c1.mul(that), Fp6Parameters);
    }

    public Fp3 mulByNonResidue(final Fp3 that) {
        return that.construct(Fp6Parameters.nonresidue().mul(that.c2), that.c0, that.c1);
    }

    public Fp6_2Over3 mul(final Fp6_2Over3 that) {
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on AbstractPairing-Friendly
     Fields.pdf; Section 3 (Karatsuba) */
        final Fp3 c0C0 = c0.mul(that.c0);
        final Fp3 c1C1 = c1.mul(that.c1);
        return new Fp6_2Over3(
                c0C0.add(mulByNonResidue(c1C1)),
                (c0.add(c1)).mul(that.c0.add(that.c1)).sub(c0C0).sub(c1C1),
                Fp6Parameters);
    }

    public Fp6_2Over3 zero() {
        return Fp6Parameters.ZERO();
    }

    public boolean isZero() {
        return c0.isZero() && c1.isZero();
    }

    public Fp6_2Over3 one() {
        return Fp6Parameters.ONE();
    }

    public boolean isOne() {
        return c0.isOne() && c1.isZero();
    }

    public Fp6_2Over3 random(final Long seed, final byte[] secureSeed) {
        return new Fp6_2Over3(c0.random(seed, secureSeed), c1.random(seed, secureSeed), Fp6Parameters);
    }

    public Fp6_2Over3 negate() {
        return new Fp6_2Over3(c0.negate(), c1.negate(), Fp6Parameters);
    }

    public Fp6_2Over3 square() {
    /* Devegili OhEig Scott Dahab --- Multiplication and Squaring on AbstractPairing-Friendly
     Fields.pdf; Section 3 (Complex squaring) */
        final Fp3 c0c1 = c0.mul(c1);
        final Fp3 factor = (c0.add(c1)).mul(c0.add(mulByNonResidue(c1)));
        return new Fp6_2Over3(
                factor.sub(c0c1).sub(mulByNonResidue(c0c1)),
                c0c1.add(c0c1),
                Fp6Parameters);
    }

    public Fp6_2Over3 inverse() {
    /* From "High-Speed Software Implementation of the Optimal Ate AbstractPairing over
    Barreto-Naehrig Curves"; Algorithm 8 */
        final Fp3 t0 = c0.square();
        final Fp3 t1 = c1.square();
        final Fp3 t2 = t0.sub(mulByNonResidue(t1));
        final Fp3 t3 = t2.inverse();

        return new Fp6_2Over3(c0.mul(t3), c1.mul(t3).negate(), Fp6Parameters);
    }

    public Fp6_2Over3 FrobeniusMap(long power) {
        return new Fp6_2Over3(
                c0.FrobeniusMap(power),
                c1.FrobeniusMap(power).mul(Fp6Parameters.FrobeniusMapCoefficientsC1()[(int) (power % 6)]),
                Fp6Parameters);
    }

    public int bitSize() {
        return Math.max(c0.bitSize(), c1.bitSize());
    }

    public Fp6_2Over3 construct(final Fp3 c0, final Fp3 c1) {
        return new Fp6_2Over3(c0, c1, Fp6Parameters);
    }

    public String toString() {
        return c0.toString() + " / " + c1.toString();
    }

    public boolean equals(final Fp6_2Over3 that) {
        if (that == null) {
            return false;
        }

        return c0.equals(that.c0) && c1.equals(that.c1);
    }
}
