/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package zk_proof_systems.zkSNARK;

import algebra.curves.barreto_naehrig.*;
import algebra.curves.barreto_naehrig.abstract_bn_parameters.AbstractBNG1Parameters;
import algebra.curves.barreto_naehrig.abstract_bn_parameters.AbstractBNG2Parameters;
import algebra.curves.barreto_naehrig.abstract_bn_parameters.AbstractBNGTParameters;
import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFr;
import algebra.curves.barreto_naehrig.bn254a.BN254aG1;
import algebra.curves.barreto_naehrig.bn254a.BN254aG2;
import algebra.curves.barreto_naehrig.bn254a.BN254aPairing;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG1Parameters;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG2Parameters;
import algebra.curves.barreto_naehrig.bn254b.BN254bFields;
import algebra.curves.barreto_naehrig.bn254b.BN254bG1;
import algebra.curves.barreto_naehrig.bn254b.BN254bG2;
import algebra.curves.barreto_naehrig.bn254b.BN254bPairing;
import algebra.curves.barreto_naehrig.bn254b.bn254b_parameters.BN254bG1Parameters;
import algebra.curves.barreto_naehrig.bn254b.bn254b_parameters.BN254bG2Parameters;
import algebra.curves.fake.*;
import algebra.curves.fake.fake_parameters.FakeFqParameters;
import algebra.curves.fake.fake_parameters.FakeG1Parameters;
import algebra.curves.fake.fake_parameters.FakeG2Parameters;
import algebra.fields.Fp;
import configuration.Configuration;
import org.junit.Before;
import org.junit.Test;
import profiler.generation.R1CSConstruction;
import relations.objects.Assignment;
import relations.r1cs.R1CSRelation;
import scala.Tuple3;
import zk_proof_systems.zkSNARK.objects.CRS;
import zk_proof_systems.zkSNARK.objects.Proof;

import java.io.Serializable;

import static org.junit.Assert.assertTrue;

public class SerialzkSNARKTest implements Serializable {
    private Configuration config;

    @Before
    public void setUp() {
        config = new Configuration();
        config.setRuntimeFlag(false);
        config.setDebugFlag(false);
    }

    private <BNFrT extends BNFields.BNFr<BNFrT>,
            BNFqT extends BNFields.BNFq<BNFqT>,
            BNFq2T extends BNFields.BNFq2<BNFqT, BNFq2T>,
            BNFq6T extends BNFields.BNFq6<BNFqT, BNFq2T, BNFq6T>,
            BNFq12T extends BNFields.BNFq12<BNFqT, BNFq2T, BNFq6T, BNFq12T>,
            BNG1T extends BNG1<BNFrT, BNFqT, BNG1T, BNG1ParametersT>,
            BNG2T extends BNG2<BNFrT, BNFqT, BNFq2T, BNG2T, BNG2ParametersT>,
            BNGTT extends BNGT<BNFqT, BNFq2T, BNFq6T, BNFq12T, BNGTT, BNGTParametersT>,
            BNG1ParametersT extends AbstractBNG1Parameters<BNFrT, BNFqT, BNG1T, BNG1ParametersT>,
            BNG2ParametersT extends AbstractBNG2Parameters<BNFrT, BNFqT, BNFq2T, BNG2T, BNG2ParametersT>,
            BNGTParametersT extends AbstractBNGTParameters<BNFqT, BNFq2T, BNFq6T, BNFq12T, BNGTT, BNGTParametersT>,
            BNPublicParametersT extends BNPublicParameters<BNFqT, BNFq2T, BNFq6T, BNFq12T>,
            BNPairingT extends BNPairing<BNFrT, BNFqT, BNFq2T, BNFq6T, BNFq12T, BNG1T, BNG2T, BNGTT, BNG1ParametersT, BNG2ParametersT, BNGTParametersT, BNPublicParametersT>>
    void SerialBNProofSystemTest(
            final int numInputs,
            final int numConstraints,
            BNFrT fieldFactory,
            BNG1T g1Factory,
            BNG2T g2Factory,
            BNPairingT pairing) {
        final Tuple3<R1CSRelation<BNFrT>, Assignment<BNFrT>, Assignment<BNFrT>> construction =
                R1CSConstruction.serialConstruct(numConstraints, numInputs, fieldFactory, config);
        final R1CSRelation<BNFrT> r1cs = construction._1();
        final Assignment<BNFrT> primary = construction._2();
        final Assignment<BNFrT> fullAssignment = construction._3();

        final CRS<BNFrT, BNG1T, BNG2T, BNGTT> CRS =
                SerialSetup.generate(r1cs, fieldFactory, g1Factory, g2Factory, pairing, config);


        final Proof<BNG1T, BNG2T> proof =
                SerialProver.prove(CRS.provingKey(), primary, fullAssignment, fieldFactory, config);

        // final boolean isValid = Verifier.verify(CRS.verificationKey(), primary, proof, pairing, config);

        // System.out.println(isValid);
        // assertTrue(isValid);
    }

    @Test
    public void SerialFakeProofSystemTest() {
        final int numInputs = 123;
        final int numConstraints = 1000000;

        FakeInitialize.init();
        final Fp fieldFactory = new FakeFqParameters().ONE();
        final FakeG1 fakeG1Factory = new FakeG1Parameters().ONE();
        final FakeG2 fakeG2Factory = new FakeG2Parameters().ONE();
        final FakePairing fakePairing = new FakePairing();

        final Tuple3<R1CSRelation<Fp>, Assignment<Fp>, Assignment<Fp>> construction = R1CSConstruction
                .serialConstruct(numConstraints, numInputs, fieldFactory, config);
        final R1CSRelation<Fp> r1cs = construction._1();
        final Assignment<Fp> primary = construction._2();
        final Assignment<Fp> auxiliary = construction._3();

        final CRS<Fp, FakeG1, FakeG2, FakeGT> CRS = SerialSetup
                .generate(r1cs, fieldFactory, fakeG1Factory, fakeG2Factory, fakePairing, config);
        final Proof<FakeG1, FakeG2> proof = SerialProver
                .prove(CRS.provingKey(), primary, auxiliary, fieldFactory, config);
        final boolean isValid = Verifier
                .verify(CRS.verificationKey(), primary, proof, fakePairing, config);

        System.out.println(isValid);
        assertTrue(isValid);
    }

//     @Test
//     public void SerialBN254aProofSystemTest() {
//         final int numInputs = 1023;
//         final int numConstraints = 1024;
//         final BN254aFr fieldFactory = BN254aFr.ONE;
//         final BN254aG1 g1Factory = BN254aG1Parameters.ONE;
//         final BN254aG2 g2Factory = BN254aG2Parameters.ONE;
//         final BN254aPairing pairing = new BN254aPairing();

//         SerialBNProofSystemTest(numInputs, numConstraints, fieldFactory, g1Factory, g2Factory, pairing);
//     }

//     @Test
//     public void SerialBN254bProofSystemTest() {
//         final int numInputs = 1023;
//         final int numConstraints = 1024;
//         final BN254bFields.BN254bFr fieldFactory = new BN254bFields.BN254bFr(1);
//         final BN254bG1 g1Factory = BN254bG1Parameters.ONE;
//         final BN254bG2 g2Factory = BN254bG2Parameters.ONE;
//         final BN254bPairing pairing = new BN254bPairing();

//         SerialBNProofSystemTest(numInputs, numConstraints, fieldFactory, g1Factory, g2Factory, pairing);
//     }
}

