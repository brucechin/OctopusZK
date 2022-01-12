/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.curves.fake.fake_parameters;

import algebra.curves.fake.FakeGT;
import algebra.curves.fake.abstract_fake_parameters.AbstractFakeGTParameters;
import algebra.math.BigInteger;

import java.io.Serializable;

public class FakeGTParameters extends AbstractFakeGTParameters implements Serializable {

    private FakeFqParameters FqParameters;

    private FakeGT ONE;

    public FakeFqParameters FqParameters() {
        if (FqParameters == null) {
            FqParameters = new FakeFqParameters();
        }

        return FqParameters;
    }

    public FakeGT ONE() {
        if (ONE == null) {
            ONE = new FakeGT(BigInteger.ONE, this);
        }

        return ONE;
    }
}
