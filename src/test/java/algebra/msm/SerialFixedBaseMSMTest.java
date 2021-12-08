/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.msm;

import algebra.fields.Fp;
import algebra.fields.fieldparameters.LargeFpParameters;
import algebra.groups.AdditiveIntegerGroup;
import algebra.groups.integergroupparameters.LargeAdditiveIntegerGroupParameters;
import org.junit.Before;
import org.junit.Test;

import java.io.Serializable;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.assertTrue;

public class SerialFixedBaseMSMTest implements Serializable {
    private LargeAdditiveIntegerGroupParameters GroupParameters;

    @Before
    public void setUp() {
        GroupParameters = new LargeAdditiveIntegerGroupParameters();
    }

    @Test
    public void NaiveMSMTest() {
        final AdditiveIntegerGroup base = new AdditiveIntegerGroup(7, GroupParameters);
        ArrayList<BigInteger> scalars = new ArrayList<>(4);
        scalars.add(new BigInteger("3"));
        scalars.add(new BigInteger("11"));
        scalars.add(new BigInteger("2"));
        scalars.add(new BigInteger("8"));

        ArrayList<AdditiveIntegerGroup> result = NaiveMSM.fixedBaseMSM(scalars, base);
        ArrayList<AdditiveIntegerGroup> answers = new ArrayList<>(4);
        answers.add(new AdditiveIntegerGroup(21, GroupParameters));
        answers.add(new AdditiveIntegerGroup(77, GroupParameters));
        answers.add(new AdditiveIntegerGroup(14, GroupParameters));
        answers.add(new AdditiveIntegerGroup(56, GroupParameters));

        for (int i = 0; i < answers.size(); i++) {
            System.out.println(result.get(i).toString() + " == " + answers.get(i).toString());
            assertTrue(result.get(i).equals(answers.get(i)));
        }
    }

    @Test
    public void SerialMSMTest() {
        final AdditiveIntegerGroup base = new AdditiveIntegerGroup(7, GroupParameters);
        final LargeFpParameters FpParameters = new LargeFpParameters();
        final Fp scalar = new Fp("4294967296", FpParameters);

        final int scalarSize = scalar.bitSize();
        final int windowSize = 2;

        List<List<AdditiveIntegerGroup>> windowTable = FixedBaseMSM
                .getWindowTable(base, scalarSize, windowSize);
        AdditiveIntegerGroup result = FixedBaseMSM
                .serialMSM(scalarSize, windowSize, windowTable, scalar);
        AdditiveIntegerGroup answers = new AdditiveIntegerGroup(1400, GroupParameters);

        System.out.println(result.toString() + " == " + answers.toString());
        assertTrue(result.equals(answers));
    }
}