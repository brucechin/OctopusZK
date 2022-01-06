package algebra.fft;

import algebra.fields.AbstractFieldElementExpanded;
import common.Combiner;
import common.MathUtils;
import common.Utils;
import configuration.Configuration;
import org.apache.spark.api.java.JavaPairRDD;
import scala.Tuple2;
import java.math.BigInteger;
import java.util.Arrays;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import org.apache.spark.TaskContext;

public class FFTAuxiliary {


    public static String byteToString(byte[] b) {
        byte[] masks = { -128, 64, 32, 16, 8, 4, 2, 1 };
        StringBuilder builder = new StringBuilder();
        builder.append('|');
        for( byte bb : b){
            for (byte m : masks) {
                if ((bb & m) == m) {
                    builder.append('1');
                } else {
                    builder.append('0');
                }
            }
            builder.append('|');
        }

        return builder.toString();
    }

    public static byte[] bigIntegerToByteArrayHelperCGBN(BigInteger bigint){
        byte[] temp = bigint.toByteArray();

        byte[] res = new byte[(temp.length + 3)/ 4 * 4];
        int new_len = (temp.length + 3) / 4 * 4;
        for(int i = 0; i < temp.length; i++){
            res[temp.length - i - 1] = temp[i];
        }
        return res;

    }

    static native <FieldT extends AbstractFieldElementExpanded<FieldT>> byte[] serialRadix2FFTNativeHelper(
        final List<byte[]> input,
        final byte[] omega, final int taskID);
    /**
     * Compute the radix-2 FFT of the vector a over the set S={omega^{0},...,omega^{m-1}}. Result is
     * stored in ArrayList input.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> void serialRadix2FFT(
            final List<FieldT> input,
            final FieldT omega, final int taskID) {
        final int n = input.size();
        final int logn = MathUtils.log2(n);
        if (n == 1) {
            return;
        }
        //System.out.println("on java side, serialRadix2FFT omega and input class type is "+ omega.getClass().getName());

        assert (n == (1 << logn));

        //--------------------------------------Java Native code--------------------------------

        ArrayList<byte[]> inputByteArray = new ArrayList<byte[]>();
        for(FieldT f : input){
            //System.out.println("JNI before=" + byteToString(f.toBigInteger().toByteArray()));
            inputByteArray.add(bigIntegerToByteArrayHelperCGBN(f.toBigInteger()));
        }

        //System.out.println("on java side serialRadix2FFT omega=" + byteToString(omega.toBigInteger().toByteArray()));
        byte[] omegaBytes = bigIntegerToByteArrayHelperCGBN(omega.toBigInteger());
        byte[] resultByteArray = FFTAuxiliary.serialRadix2FFTNativeHelper(inputByteArray, omegaBytes, taskID);

        int size_of_bigint_cpp_side = 64;
        for(int i = 0; i < input.size(); i++){
            byte[] slice = Arrays.copyOfRange(resultByteArray, i*size_of_bigint_cpp_side, (i+1)*size_of_bigint_cpp_side);
            byte[] converted_back = new byte[64];
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back[j] = slice[size_of_bigint_cpp_side - j - 1];
            }
            BigInteger bi = new BigInteger(converted_back);
            //System.out.println("JNI after =" + byteToString(bi.toByteArray()));
            FieldT temp = input.get(0).zero();
            temp.setBigInteger(bi);
            input.set(i, temp);

        }

        //--------------------------------------Java Native code--------------------------------

    

    }

    /**
     * A distributed version of serialRadix2FFT.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, FieldT>
    distributedRadix2FFT(
            final JavaPairRDD<Long, FieldT> input,
            final long rows,
            final long columns,
            final boolean inverse,
            final FieldT fieldFactory) {
        assert (MathUtils.isPowerOfTwo(rows));
        assert (MathUtils.isPowerOfTwo(columns));

        /* Initialization */
        final long size = rows * columns;
        final FieldT omegaShift = fieldFactory.rootOfUnity(size);
        final Combiner<FieldT> combine = new Combiner<>();
        final SerialFFT<FieldT> rowDomain = new SerialFFT<>(rows, fieldFactory);
        final SerialFFT<FieldT> columnDomain = new SerialFFT<>(columns, fieldFactory);

        /* Algorithm 1: Forward FFT, Mapper */
        final JavaPairRDD<Long, FieldT> columnGroups = input.mapToPair(element -> {
            /* AbstractGroup the array of inputs into rows using the combiner. */
            final long group = element._1 % rows;
            final long index = element._1 / rows;

            return new Tuple2<>(group, new Tuple2<>(index, element._2));
        }).combineByKey(combine.createGroup, combine.mergeElement, combine.mergeCombiner)
                .mapValues(partition -> {
                    TaskContext tc = TaskContext.get();
                    long taskID = tc.taskAttemptId();
                    /* Compute the forward FFT, on domain size columns, for each group. */
                    ArrayList<FieldT> groupArray = Utils.convertFromPairs(partition, (int) columns);

                    if (inverse) {
                        columnDomain.radix2InverseFFT(groupArray, (int)taskID);
                    } else {
                        columnDomain.radix2FFT(groupArray, (int)taskID);
                    }

                    return groupArray;
                }).flatMapToPair(element -> {
                    /* Bitshift and map to key j. */
                    final long index = element._1;

                    ArrayList<Tuple2<Long, FieldT>> combinedNumbers = new ArrayList<>();
                    for (int i = 0; i < columns; i++) {
                        final FieldT nthRoot =
                                inverse ? omegaShift.pow(index * i).inverse() : omegaShift.pow(index * i);
                        combinedNumbers.add(new Tuple2<>(i * rows + index, nthRoot.mul(element._2.get(i))));
                    }

                    return combinedNumbers.iterator();
                });

        /* Algorithm 2: Forward FFT, Reducer */
        return columnGroups.mapToPair(element -> {
            /* AbstractGroup the array of inputs into columns using the combiner. */
            final long group = element._1 / rows;
            final long index = element._1 % rows;

            return new Tuple2<>(group, new Tuple2<>(index, element._2));
        }).combineByKey(combine.createGroup, combine.mergeElement, combine.mergeCombiner)
                .mapValues(partition -> {
                    TaskContext tc = TaskContext.get();
                    long taskID = tc.taskAttemptId();
                    /* Compute the forward FFT, on domain size rows, for each group. */
                    ArrayList<FieldT> groupArray = Utils.convertFromPairs(partition, (int) rows);

                    if (inverse) {
                        rowDomain.radix2InverseFFT(groupArray, (int)taskID);
                    } else {
                        rowDomain.radix2FFT(groupArray, (int)taskID);
                    }

                    return groupArray;
                }).flatMapToPair(element -> {
                    /* Serialize and order evaluation results. */
                    final long index = element._1;

                    ArrayList<Tuple2<Long, FieldT>> outputs = new ArrayList<>();
                    for (int i = 0; i < rows; i++) {
                        outputs.add(new Tuple2<>(i * columns + index, element._2.get(i)));
                    }

                    return outputs.iterator();
                });
    }

    /**
     * Translate the vector input to a coset defined by g. Result is stored in ArrayList input.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> void multiplyByCoset(
            final List<FieldT> input,
            final FieldT g) {
        FieldT coset = g;
        for (int i = 1; i < input.size(); ++i) {
            input.set(i, input.get(i).mul(coset));
            coset = coset.mul(g);
        }
    }

    /**
     * A distributed version of multiplyByCoset.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, FieldT>
    distributedMultiplyByCoset(
            final JavaPairRDD<Long, FieldT> input,
            final FieldT g) {
        //TODO lianke implement this function in GPU
        return input.mapToPair(term -> new Tuple2<>(term._1, term._2.mul(g.pow(term._1))));
    }

    /**
     * Compute the m Lagrange coefficients, relative to the set S={omega^{0},...,omega^{m-1}}, at the
     * field element t.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> List<FieldT>
    serialRadix2LagrangeCoefficients(
            final FieldT t,
            final int m) {
        if (m == 1) {
            final List<FieldT> coeff = new ArrayList<>();
            coeff.add(t.one());
            return coeff;
        }

        assert (m == (1 << MathUtils.log2(m)));

        final FieldT omega = t.rootOfUnity(m);
        final FieldT one = t.one();

        List<FieldT> lagrangeCoefficients = new ArrayList<>(Collections.nCopies(m, t.zero()));

        /*
         * If t equals one of the roots of unity in S={omega^{0},...,omega^{m-1}}
         * then output 1 at the right place, and 0 elsewhere
         */
        if (t.pow(m).equals(one)) {
            FieldT omega_i = one;
            for (int i = 0; i < m; i++) {
                /* i.e., t equals omega^i */
                if (omega_i.equals(t)) {
                    lagrangeCoefficients.set(i, one);
                    return lagrangeCoefficients;
                }

                omega_i = omega_i.mul(omega);
            }
        }

        /*
         * Otherwise, if t does not equal any of the roots of unity in S,
         * then compute each L_{i,S}(t) as Z_{S}(t) * v_i / (t-\omega^i)
         * where:
         * - Z_{S}(t) = \prod_{j} (t-\omega^j) = (t^m-1), and
         * - v_{i} = 1 / \prod_{j \neq i} (\omega^i-\omega^j).
         * Below we use the fact that v_{0} = 1/m and v_{i+1} = \omega * v_{i}.
         */
        final FieldT Z = t.pow(m).sub(one);
        final FieldT mInverse = t.construct(m).inverse();
        FieldT l = Z.mul(mInverse);
        FieldT r = one;
        for (int i = 0; i < m; ++i) {
            lagrangeCoefficients.set(i, l.mul(t.sub(r).inverse()));
            l = l.mul(omega);
            r = r.mul(omega);
        }

        return lagrangeCoefficients;
    }

    /**
     * A distributed version of serialRadix2LagrangeCoefficients.
     */
    static <FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, FieldT>
    distributedRadix2LagrangeCoefficients(
            final FieldT t,
            final long m,
            final Configuration config) {
        if (m == 1) {
            return config.sparkContext()
                    .parallelizePairs(Collections.singletonList(new Tuple2<>((long) 0, t.one())));
        }

        assert (m == (1 << MathUtils.log2(m)));

        final FieldT omega = t.rootOfUnity(m);
        final FieldT one = t.one();
        final FieldT zero = t.zero();

        /*
         * If t equals one of the roots of unity in S={omega^{0},...,omega^{m-1}}
         * then output 1 at the right place, and 0 elsewhere
         */
        if (t.pow(m).equals(one)) {
            return Utils.fillRDD(m, zero, config).mapToPair(term -> {
                if (!omega.pow(term._1).equals(t)) {
                    return term;
                }
                return new Tuple2<>(term._1, one);
            });
        }

        /*
         * Otherwise, if t does not equal any of the roots of unity in S,
         * then compute each L_{i,S}(t) as Z_{S}(t) * v_i / (t-\omega^i)
         * where:
         * - Z_{S}(t) = \prod_{j} (t-\omega^j) = (t^m-1), and
         * - v_{i} = 1 / \prod_{j \neq i} (\omega^i-\omega^j).
         * Below we use the fact that v_{0} = 1/m and v_{i+1} = \omega * v_{i}.
         */
        final FieldT Z = t.pow(m).sub(one);
        final FieldT mInverse = t.construct(m).inverse();
        final FieldT l = Z.mul(mInverse);

        return Utils.fillRDD(m, zero, config).mapToPair(term -> {
            final FieldT omega_i = omega.pow(term._1);
            final FieldT l_i = l.mul(omega_i);

            return new Tuple2<>(term._1, l_i.mul(t.sub(omega_i).inverse()));
        });
    }

    public static <FieldT extends AbstractFieldElementExpanded<FieldT>> HashMap<Long, FieldT>
    subsequenceRadix2LagrangeCoefficients(
            final FieldT t,
            final long m,
            final List<Long> indices) {

        assert (m == (1 << MathUtils.log2(m)));

        final HashMap<Long, FieldT> coefficients = new HashMap<>();

        final FieldT omega = t.rootOfUnity(m);
        final FieldT one = t.one();
        final FieldT zero = t.zero();

        /*
         * If t equals one of the roots of unity in S={omega^{0},...,omega^{m-1}}
         * then output 1 at the right place, and 0 elsewhere
         */
        if (t.pow(m).equals(one)) {
            for (Long index : indices) {
                if (omega.pow(index).equals(t)) {
                    coefficients.put(index, one);
                } else {
                    coefficients.put(index, zero);
                }
            }
            return coefficients;
        }

        /*
         * Otherwise, if t does not equal any of the roots of unity in S,
         * then compute each L_{i,S}(t) as Z_{S}(t) * v_i / (t-\omega^i)
         * where:
         * - Z_{S}(t) = \prod_{j} (t-\omega^j) = (t^m-1), and
         * - v_{i} = 1 / \prod_{j \neq i} (\omega^i-\omega^j).
         * Below we use the fact that v_{0} = 1/m and v_{i+1} = \omega * v_{i}.
         */
        final FieldT Z = t.pow(m).sub(one);
        final FieldT mInverse = t.construct(m).inverse();
        final FieldT l = Z.mul(mInverse);

        for (Long index : indices) {
            final FieldT omega_i = omega.pow(index);
            final FieldT l_i = l.mul(omega_i);

            coefficients.put(index, l_i.mul(t.sub(omega_i).inverse()));
        }
        return coefficients;
    }
}
