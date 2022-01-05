package profiler.profiling;

import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFr;
import algebra.curves.barreto_naehrig.bn254a.BN254aG1;
import algebra.curves.barreto_naehrig.bn254a.BN254aG2;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG1Parameters;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG2Parameters;
import algebra.msm.VariableBaseMSM;
import configuration.Configuration;
import org.apache.spark.api.java.JavaRDD;
import profiler.generation.VariableBaseMSMGenerator;
import scala.Tuple2;

import java.util.ArrayList;

public class VariableBaseMSMProfiling {

    public static void serialVariableBaseMSMG1Profiling(final Configuration config, final long size) throws Exception{
        final BN254aG1 groupFactory = new BN254aG1Parameters().ONE();
        final BN254aFr fieldFactory = new BN254aFr(2L);

        final ArrayList<BN254aFr> scalars = new ArrayList<>();
        final ArrayList<BN254aG1> bases = new ArrayList<>();

        for (long j = 0; j < size; j++) {
            scalars.add(fieldFactory.random(config.seed(), config.secureSeed()));
            bases.add(groupFactory);
        }
        System.out.println("random benchmark data generation done.");

        config.setContext("VariableBaseMSMG1-Serial");
        config.beginRuntimeMetadata("Size (inputs)", size);

        config.beginLog("VariableBaseMSM");
        config.beginRuntime("VariableBaseMSM");
        final BN254aG1 res = VariableBaseMSM.serialMSM(scalars, bases);
        config.endRuntime("VariableBaseMSM");
        config.endLog("VariableBaseMSM");

        config.writeRuntimeLog(config.context());
        System.out.println("VarMSM output =" +res.toString());

    }

    public static void serialVariableBaseMSMG2Profiling(final Configuration config, final long size) throws Exception{
        final BN254aG2 groupFactory = new BN254aG2Parameters().ONE();
        final BN254aFr fieldFactory = new BN254aFr(2L);

        final ArrayList<BN254aFr> scalars = new ArrayList<>();
        final ArrayList<BN254aG2> bases = new ArrayList<>();

        for (long j = 0; j < size; j++) {
            scalars.add(fieldFactory.random(config.seed(), config.secureSeed()));
            bases.add(groupFactory);
        }
        System.out.println("random benchmark data generation done.");

        config.setContext("VariableBaseMSMG2-Serial");
        config.beginRuntimeMetadata("Size (inputs)", size);

        config.beginLog("VariableBaseMSM");
        config.beginRuntime("VariableBaseMSM");
        final BN254aG2 res = VariableBaseMSM.serialMSM(scalars, bases);
        System.out.println("res="+res.toString());
        config.endRuntime("VariableBaseMSM");
        config.endLog("VariableBaseMSM");
        System.out.println("VarMSM output =" +res.toString());
        config.writeRuntimeLog(config.context());
    }

    public static void distributedVariableBaseMSMG1Profiling(
            final Configuration config,
            final long size) throws Exception{
        final JavaRDD<Tuple2<BN254aFr, BN254aG1>> input = VariableBaseMSMGenerator
                .generateG1Data(config, size);

        config.setContext("VariableBaseMSMG1");
        config.beginRuntimeMetadata("Size (inputs)", size);
        System.out.println("random benchmark data generation done.");

        config.beginLog("VariableBaseMSM");
        config.beginRuntime("VariableBaseMSM");
        final BN254aG1 res2 = VariableBaseMSM.distributedMSM(input);
        System.out.println("res="+res2.toString());
        config.endRuntime("VariableBaseMSM");
        config.endLog("VariableBaseMSM");

        config.writeRuntimeLog(config.context());
    }

    public static void distributedVariableBaseMSMG2Profiling(
            final Configuration config,
            final long size) throws Exception{
        final JavaRDD<Tuple2<BN254aFr, BN254aG2>> input = VariableBaseMSMGenerator
                .generateG2Data(config, size);

        config.setContext("VariableBaseMSMG2");
        config.beginRuntimeMetadata("Size (inputs)", size);
        System.out.println("random benchmark data generation done.");

        config.beginLog("VariableBaseMSM");
        config.beginRuntime("VariableBaseMSM");
        final BN254aG2 res2 = VariableBaseMSM.distributedMSM(input);
        System.out.println("res="+res2.toString());
        config.endRuntime("VariableBaseMSM");
        config.endLog("VariableBaseMSM");

        config.writeRuntimeLog(config.context());
    }


}
