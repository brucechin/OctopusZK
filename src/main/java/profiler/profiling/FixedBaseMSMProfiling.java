package profiler.profiling;

import algebra.curves.barreto_naehrig.bn254a.BN254aFields.BN254aFr;
import algebra.curves.barreto_naehrig.bn254a.BN254aG1;
import algebra.curves.barreto_naehrig.bn254a.BN254aG2;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG1Parameters;
import algebra.curves.barreto_naehrig.bn254a.bn254a_parameters.BN254aG2Parameters;
import algebra.msm.FixedBaseMSM;
import configuration.Configuration;
import org.apache.spark.api.java.JavaPairRDD;
import profiler.generation.FixedBaseMSMGenerator;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

public class FixedBaseMSMProfiling {

    public static void serialFixedBaseMSMG1Profiling(final Configuration config, final long size) throws Exception {
        final BN254aFr fieldFactory = new BN254aFr(2L);
        final BN254aG1 groupFactory = new BN254aG1Parameters().ONE();
        BN254aG1 generatorG1 = groupFactory.random(config.seed(), config.secureSeed());
        final int scalarSizeG1 = generatorG1.bitSize();

        final Random rand = new Random(System.nanoTime());
        final ArrayList<BN254aFr> scalars = new ArrayList<>();
        for (long i = 0; i < size; i++) {
            scalars.add(fieldFactory.random(rand.nextLong(), null));
        }
        System.out.println("random benchmark data generation done.");

        config.setContext("FixedBaseMSMG1-Serial");
        config.beginRuntimeMetadata("Size (inputs)", size);

        final int windowSizeG1 = FixedBaseMSM.getWindowSize((int) size, generatorG1);

        final int numWindowsG1 = (scalarSizeG1 % windowSizeG1 == 0) ? scalarSizeG1 / windowSizeG1 : scalarSizeG1 / windowSizeG1+ 1;
        final int innerLimitG1 = (int) Math.pow(2, windowSizeG1);

        config.beginLog("FixedBaseMSM");
        config.beginRuntime("FixedBaseMSM");
        final List<BN254aG1> result = FixedBaseMSM
                .batchMSM(scalarSizeG1, windowSizeG1, numWindowsG1, innerLimitG1, generatorG1, scalars);
        config.endRuntime("FixedBaseMSM");
        config.endLog("FixedBaseMSM");

        config.writeRuntimeLog(config.context());
        System.out.println("FixMSM output =" +result.get(result.size() - 1).toString());
    }

    public static void serialFixedBaseMSMG2Profiling(final Configuration config, final long size) throws Exception{
        final BN254aFr fieldFactory = new BN254aFr(2L);
        BN254aG2 groupFactory = new BN254aG2Parameters().ONE();
        BN254aG2 generatorG2 = groupFactory.random(config.seed(), config.secureSeed());

        final int scalarSizeG2 = generatorG2.bitSize();

        final Random rand = new Random(System.nanoTime());
        final ArrayList<BN254aFr> scalars = new ArrayList<>();
        for (long i = 0; i < size; i++) {
            scalars.add(fieldFactory.random(rand.nextLong(), null));
        }
        System.out.println("random benchmark data generation done.");
        config.setContext("FixedBaseMSMG2-Serial");
        config.beginRuntimeMetadata("Size (inputs)", size);

        config.beginLog("FixedBaseMSM");
        config.beginRuntime("FixedBaseMSM");
        final int windowSizeG2 = FixedBaseMSM.getWindowSize((int) size, generatorG2);
        final int numWindowsG2 = (scalarSizeG2 % windowSizeG2 == 0) ? scalarSizeG2 / windowSizeG2 : scalarSizeG2 / windowSizeG2+ 1;
        final int innerLimitG2 = (int) Math.pow(2, windowSizeG2);


        final List<BN254aG2> result = FixedBaseMSM
        .batchMSM(scalarSizeG2, windowSizeG2, numWindowsG2, innerLimitG2, generatorG2, scalars);
        
        config.endRuntime("FixedBaseMSM");
        config.endLog("FixedBaseMSM");

        config.writeRuntimeLog(config.context());
        System.out.println("FixMSM output =" +result.get(result.size() - 1).toString());

    }

    public static void distributedFixedBaseMSMG1Profiling(final Configuration config, final long size) {
        final BN254aG1 groupFactory = new BN254aG1Parameters().ONE();
        BN254aG1 generatorG1 = groupFactory.random(config.seed(), config.secureSeed());
        final int scalarSizeG1 = generatorG1.bitSize();

        final JavaPairRDD<Long, BN254aFr> scalars = FixedBaseMSMGenerator.generateData(config, size);
        System.out.println("random benchmark data generation done.");

        config.setContext("FixedBaseMSMG1");
        config.beginRuntimeMetadata("Size (inputs)", size);
        config.beginLog("Finished random R1CS generation for benchmark purpose. benchmark data generation is not counted into the setup and proving time.");


        final int windowSizeG1 = FixedBaseMSM
                .getWindowSize(size / config.numPartitions(), generatorG1);
        final int numWindowsG1 = (scalarSizeG1 % windowSizeG1 == 0) ? scalarSizeG1 / windowSizeG1 : scalarSizeG1 / windowSizeG1+ 1;
        final int innerLimitG1 = (int) Math.pow(2, windowSizeG1);
        config.beginLog("FixedBaseMSM");
        config.beginRuntime("FixedBaseMSM");
        long start = System.currentTimeMillis();
        JavaPairRDD<Long, BN254aG1> result = FixedBaseMSM.distributedBatchMSM(
                scalarSizeG1,
                windowSizeG1,
                numWindowsG1, innerLimitG1,
                generatorG1,
                scalars,
                config.sparkContext()).persist(config.storageLevel());
        result.count();
        long finish = System.currentTimeMillis();
        long timeElapsed = finish - start;
        System.out.println("distributedFixedBaseMSMG1Profiling size="+ size +  "  time elapsed: " + timeElapsed + " ms");

        config.endRuntime("FixedBaseMSM");
        config.endLog("FixedBaseMSM");

        config.writeRuntimeLog(config.context());
    }

    public static void distributedFixedBaseMSMG2Profiling(final Configuration config, final long size) {

        BN254aG2 groupFactory = new BN254aG2Parameters().ONE();
        BN254aG2 generatorG2 = groupFactory.random(config.seed(), config.secureSeed());
        final int scalarSizeG2 = generatorG2.bitSize();

        final JavaPairRDD<Long, BN254aFr> scalars = FixedBaseMSMGenerator.generateData(config, size);
        System.out.println("random benchmark data generation done.");

        config.setContext("FixedBaseMSMG2");
        config.beginRuntimeMetadata("Size (inputs)", size);

        config.beginLog("FixedBaseMSM");
        config.beginRuntime("FixedBaseMSM");
        final int windowSizeG2 = FixedBaseMSM.getWindowSize((int) size, generatorG2);
        final int numWindowsG2 = (scalarSizeG2 % windowSizeG2 == 0) ? scalarSizeG2 / windowSizeG2 : scalarSizeG2 / windowSizeG2+ 1;
        final int innerLimitG2 = (int) Math.pow(2, windowSizeG2);

        FixedBaseMSM.distributedBatchMSM(
                scalarSizeG2,
                windowSizeG2,
                numWindowsG2, innerLimitG2,
                generatorG2,
                scalars,
                config.sparkContext()).count();
        config.endRuntime("FixedBaseMSM");
        config.endLog("FixedBaseMSM");

        config.writeRuntimeLog(config.context());
    }
}
