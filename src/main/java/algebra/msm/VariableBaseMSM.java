/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.msm;

import algebra.fields.AbstractFieldElementExpanded;
import algebra.groups.AbstractGroup;
import common.MathUtils;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.PriorityQueue;
import org.apache.spark.api.java.JavaRDD;
import scala.Tuple2;
import java.util.Arrays;
import java.io.*;

public class VariableBaseMSM {

    public static final BigInteger BOS_COSTER_MSM_THRESHOLD = new BigInteger("1048576");

    static {
		System.loadLibrary("AlgebraMSMVariableBaseMSM");
        System.out.println("AlgebraMSMVariableBaseMSM loaded");
	}
    
    /**
     * The algorithm starts by sorting the input in order of scalar size. For each iteration, it
     * computes the difference of the two largest scalars, and multiplies the accrued base by this
     * difference of exponents. This result is then summed and lastly returned.
     */
    public static <GroupT extends AbstractGroup<GroupT>> GroupT sortedMSM(
             final List<Tuple2<BigInteger, GroupT>> input) {

        GroupT result = input.get(0)._2.zero();
        GroupT base = input.get(0)._2.zero();
        BigInteger scalar;

        input.sort((v1, v2) -> v1._1.compareTo(v2._1));

        for (int i = input.size() - 1; i >= 0; i--) {
            scalar = i != 0 ? input.get(i)._1.subtract(input.get(i - 1)._1) : input.get(i)._1;
            base = base.add(input.get(i)._2);
            result = result.add(base.mul(scalar));
        }
        return result;
    }

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



    /**
     * The Bos-Coster algorithm works by repeatedly recursively computing
     * (e1 - e2) * b1 + e2 * (b1 + b2) where the scalars are ordered such that
     * e1 >= e2 >= ... >= en in a priority queue. The result is summed and
     * returned after all scalar-base pairs have been computed.
     *
     * BigInteger values must be positive.
     */
    public static <GroupT extends AbstractGroup<GroupT>> GroupT bosCosterMSM(
             final List<Tuple2<BigInteger, GroupT>> input) {

        final PriorityQueue<Tuple2<BigInteger, GroupT>> sortedScalarPairs =
                new PriorityQueue<>(input.size(), (v1, v2) -> v2._1.compareTo(v1._1));
        sortedScalarPairs.addAll(input);

        GroupT result = input.get(0)._2.zero();
        Tuple2<BigInteger, GroupT> e1;
        Tuple2<BigInteger, GroupT> e2;

        while ((e1 = sortedScalarPairs.poll()) != null && (e2 = sortedScalarPairs.poll()) != null) {

            if (e1._1.divide(e2._1).compareTo(BOS_COSTER_MSM_THRESHOLD) >= 0) {
                result = result.add(e1._2.mul(e1._1));
                sortedScalarPairs.add(e2);
            } else {
                final BigInteger value = e1._1.subtract(e2._1);

                if (!value.equals(BigInteger.ZERO)) {
                    sortedScalarPairs.add(new Tuple2<>(value, e1._2));
                }

                sortedScalarPairs.add(new Tuple2<>(e2._1, e1._2.add(e2._2)));
            }
        }

        while (e1 != null) {
            result = result.add(e1._2.mul(e1._1));
            e1 = sortedScalarPairs.poll();
        }

        return result;
    }

    public static byte[] bigIntegerToByteArrayHelperCGBN(BigInteger bigint){
        byte[] temp = bigint.toByteArray();
        //byte[] res = new byte[(temp.length + 3)/ 4 * 4];
        byte[] res = new byte[32]; //254-bit BN254 biginteger is enough to be saved in 32bytes.
        for(int i = 0; i < temp.length; i++){
            res[temp.length - i - 1] = temp[i];
        }
        
        return res;

    }
    public static byte[] bigIntegerToByteArrayHelper(BigInteger bigint){
        byte[] temp = bigint.toByteArray();

        byte[] res = new byte[temp.length / 4 * 4 + 4];
        int size_diff = temp.length / 4 * 4 + 4 - temp.length;
        int j = temp.length - 1;
        for(; j >= 3; j-=4){
            res[ size_diff+ j] = temp[j-3];
            res[size_diff+ j-1] = temp[j-2];
            res[ size_diff+j-2] = temp[j-1];
            res[ size_diff+j-3] = temp[j];
        }
        if(j == 2){
            res[2] = temp[j-2];
            res[1] = temp[j-1];
            res[0] = temp[j];
        }else if(j == 1){
            res[1] = temp[j-1];
            res[0] = temp[j];
        }else if(j == 0){
            res[0] = temp[j];
        }
        return res;

    }


    public static <GroupT extends AbstractGroup<GroupT>> GroupT pippengerMSM(
             final List<Tuple2<BigInteger, GroupT>> input, final int numBits) {

        final int length = input.size();
        final int log2Length = Math.max(1, MathUtils.log2(length));
        final int c = log2Length - (log2Length / 3);
        System.out.println("on pippengerMSM, input class type is "+ input.get(0)._2.getClass().getName());

        final int numBuckets = 1 << c;
        final int numGroups = (numBits + c - 1) / c;

        final GroupT zero = input.get(0)._2.zero();
        final ArrayList<GroupT> bucketsModel = new ArrayList<>(Collections.nCopies(numBuckets, zero));

        GroupT result = zero;
        //System.out.println("java side length: " + length + "log2 length: " + log2Length + "c " + c + " numGroups " + numGroups + "numBuckets " +  numBuckets);
        ArrayList<GroupT> resultArray = new ArrayList<>();
        //TOD lianke: this for loop could mean we call the CUDA kernel for numGroups times.
        for (int k = numGroups - 1; k >= 0; k--) {


            //TODO lianke: this buckets should be a global array shared by all CUDA processes.
            final ArrayList<GroupT> buckets = new ArrayList<>(bucketsModel);
            //TODO lianke: this for loop should be CUDA implemented.
            for (int i = 0; i < length; i++) {
                int id = 0;
                for (int j = 0; j < c; j++) {
                    if (input.get(i)._1.testBit(k * c + j)) {
                        id |= 1 << j;
                    }
                }

                if (id == 0) {
                    continue;
                }

                // Potentially use mixed addition here.
                buckets.set(id, buckets.get(id).add(input.get(i)._2));
            }

            GroupT runningSum = zero;

            for (int i = numBuckets - 1; i > 0; i--) {
                // Potentially use mixed addition here.
                runningSum = runningSum.add(buckets.get(i));
                result = result.add(runningSum);
            }
            resultArray.add(result);

            if (k > 0) {
                for (int i = 0; i < c; i++) {
                    result = result.twice();
                }
            }
            //System.out.println("k=" +k + " result=" + result.toString());
        }

        return result;
    }

    
 

    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] variableBaseSerialMSMNativeHelper(
            final byte[] basesXYZ,
            final byte[] scalars,
            final int batch_size);

    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT serialMSM(final List<FieldT> scalars, final List<GroupT> bases) throws Exception {
        // System.out.println("variableBaseSerialMSM info:");
        // System.out.println("variableBaseSerialMSM base size :" + bases.size() + " type:" +bases.get(0).getClass().getName());
        // System.out.println("variableBaseSerialMSM scalars size :" + scalars.size() + " type:"+scalars.get(0).getClass().getName());

        assert (bases.size() == scalars.size());


        //serialMSM is only called for  BN254G1 curve.

        ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        for (FieldT scalar : scalars) {
            bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));
        }
        byte[] bigScalarByteArray =  bigScalarStream.toByteArray();



        ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );

        for (GroupT base : bases){
            ArrayList<BigInteger> three_values = base.BN254G1ToBigInteger();
            outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
            outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
            outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));

        }
        byte[] baseByteArrayXYZ = outputStream.toByteArray();

        byte[] resArray = variableBaseSerialMSMNativeHelper(baseByteArrayXYZ, bigScalarByteArray, bases.size());
        
        int size_of_bigint_cpp_side = 64;
        
        byte[] converted_back_X = new byte[64];
        byte[] converted_back_Y = new byte[64];
        byte[] converted_back_Z = new byte[64];

        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_X[j] = resArray[size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Y[j] = resArray[2*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Z[j] = resArray[3*size_of_bigint_cpp_side - j - 1];
        }
        BigInteger bi_X = new BigInteger(converted_back_X);
        BigInteger bi_Y = new BigInteger(converted_back_Y);
        BigInteger bi_Z = new BigInteger(converted_back_Z);       
        GroupT jni_res = bases.get(0).zero();
        jni_res.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);

        //System.out.println("CUDA output=" + jni_res.toString());


        return jni_res;
    }


    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] variableBaseDoubleMSMNativeHelper(
            final byte[] bases1,
            final byte[] bases2,
            final byte[] scalars,
            final int batch_size1);

    public static <
            T1 extends AbstractGroup<T1>,
            T2 extends AbstractGroup<T2>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    Tuple2<T1, T2> doubleMSM(final List<FieldT> scalars, final List<Tuple2<T1, T2>> bases) throws Exception {
        // System.out.println("variableBaseDoubleMSM info:");
        // System.out.println("variableBaseDoubleMSM base size :" + bases.size() + " type:" +bases.get(0)._1.getClass().getName() + " " +bases.get(0)._2.getClass().getName()  );
        // System.out.println("variableBaseDoubleMSM scalars size :" + scalars.size() + " type:"+scalars.get(0).getClass().getName());

        assert (bases.size() == scalars.size());

        final int size = bases.size();
        assert (size > 0);



        ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        for (FieldT scalar : scalars) {
            bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));
        }
        byte[] bigScalarByteArray =  bigScalarStream.toByteArray();


        ByteArrayOutputStream outputStreamBase1 = new ByteArrayOutputStream( );
        ByteArrayOutputStream outputStreamBase2 = new ByteArrayOutputStream( );


        for (Tuple2<T1, T2> base : bases){
            ArrayList<BigInteger> three_values = base._1.BN254G1ToBigInteger();
            outputStreamBase1.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
            outputStreamBase1.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
            outputStreamBase1.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));

            ArrayList<BigInteger> six_values = base._2.BN254G2ToBigInteger();
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));   
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
            outputStreamBase2.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));   
        }
        byte[] baseByteArrayBase1 = outputStreamBase1.toByteArray();
        byte[] baseByteArrayBase2 = outputStreamBase2.toByteArray();

        byte[] resArray = variableBaseDoubleMSMNativeHelper(baseByteArrayBase1, baseByteArrayBase2, bigScalarByteArray, bases.size());


        int size_of_bigint_cpp_side = 64;


        byte[] converted_back_X = new byte[64];
        byte[] converted_back_Y = new byte[64];
        byte[] converted_back_Z = new byte[64];

        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_X[j] = resArray[size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Y[j] = resArray[2*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Z[j] = resArray[3*size_of_bigint_cpp_side - j - 1];
        }
        BigInteger bi_X = new BigInteger(converted_back_X);
        BigInteger bi_Y = new BigInteger(converted_back_Y);
        BigInteger bi_Z = new BigInteger(converted_back_Z);       
        T1 jni_res1 = bases.get(0)._1.zero();
        jni_res1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);




        byte[] converted_back_Xa = new byte[64];
        byte[] converted_back_Ya = new byte[64];
        byte[] converted_back_Za = new byte[64];
        byte[] converted_back_Xb = new byte[64];
        byte[] converted_back_Yb = new byte[64];
        byte[] converted_back_Zb = new byte[64];

        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Xa[j] = resArray[4*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Xb[j] = resArray[5*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Ya[j] = resArray[6*size_of_bigint_cpp_side - j - 1];
        }

        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Yb[j] = resArray[7*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Za[j] = resArray[8*size_of_bigint_cpp_side - j - 1];
        }
        for(int j =0; j < size_of_bigint_cpp_side; j++){
            converted_back_Zb[j] = resArray[9*size_of_bigint_cpp_side - j - 1];
        }

        BigInteger bi_Xa = new BigInteger(converted_back_Xa);
        BigInteger bi_Ya = new BigInteger(converted_back_Ya);
        BigInteger bi_Za = new BigInteger(converted_back_Za);    
        BigInteger bi_Xb = new BigInteger(converted_back_Xb);
        BigInteger bi_Yb = new BigInteger(converted_back_Yb);
        BigInteger bi_Zb = new BigInteger(converted_back_Zb);      
        T2 jni_res2 = bases.get(0)._2.zero();
        jni_res2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);

        return new Tuple2<>(jni_res1, jni_res2);




    }

    public static <T1 extends AbstractGroup<T1>, T2 extends AbstractGroup<T2>>
    Tuple2<T1, T2> doubleMSM(final List<Tuple2<BigInteger, Tuple2<T1, T2>>> input) {
        System.out.println("doubleMSM with only one input parameter;" );
        final int size = input.size();
        assert (size > 0);

        ArrayList<Tuple2<BigInteger, T1>> converted1 = new ArrayList<>(size);
        ArrayList<Tuple2<BigInteger, T2>> converted2 = new ArrayList<>(size);

        T1 acc1 = input.get(0)._2._1.zero();
        T2 acc2 = input.get(0)._2._2.zero();

        int numBits = 0;
        for (int i = 0; i < size; i++) {
            final Tuple2<T1, T2> value = input.get(i)._2;
            final BigInteger scalar = input.get(i)._1;
            if (scalar.equals(BigInteger.ZERO)) {
                continue;
            }

            // Mixed addition
            if (scalar.equals(BigInteger.ONE)) {
                acc1 = acc1.add(value._1);
                acc2 = acc2.add(value._2);
            } else {
                converted1.add(new Tuple2<>(scalar, value._1));
                converted2.add(new Tuple2<>(scalar, value._2));
                numBits = Math.max(numBits, scalar.bitLength());
            }
        }

        return new Tuple2<>(
                converted1.isEmpty() ? acc1 : acc1.add(VariableBaseMSM.pippengerMSM(converted1, numBits)),
                converted2.isEmpty() ? acc2 : acc2.add(VariableBaseMSM.pippengerMSM(converted2, numBits)));
    }

    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT distributedMSM(final JavaRDD<Tuple2<FieldT, GroupT>> input) {

        return input.mapPartitions(partition -> {
            final List<Tuple2<BigInteger, GroupT>> pairs = new ArrayList<>();

            int numBits = 0;
            while (partition.hasNext()) {
                final Tuple2<FieldT, GroupT> pair = partition.next();

                final BigInteger scalar = pair._1.toBigInteger();
                if (scalar.equals(BigInteger.ZERO)) {
                    continue;
                }

                pairs.add(new Tuple2<>(scalar, pair._2));
                numBits = Math.max(numBits, scalar.bitLength());
            }

            return
                    pairs.isEmpty() ?
                            Collections.emptyListIterator() :
                            Collections.singletonList(pippengerMSM(pairs, numBits)).iterator();
        }).reduce(GroupT::add);
    }

    public static <
            G1T extends AbstractGroup<G1T>,
            G2T extends AbstractGroup<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    Tuple2<G1T, G2T> distributedDoubleMSM(final JavaRDD<Tuple2<FieldT, Tuple2<G1T, G2T>>> input) {

        return input.mapPartitions(partition -> {
            final List<Tuple2<BigInteger, Tuple2<G1T, G2T>>> pairs = new ArrayList<>();
            while (partition.hasNext()) {
                Tuple2<FieldT, Tuple2<G1T, G2T>> pair = partition.next();
                final BigInteger scalar = pair._1.toBigInteger();
                if (scalar.equals(BigInteger.ZERO)) {
                    continue;
                }
                pairs.add(new Tuple2<>(scalar, pair._2));
            }
            return
                    pairs.isEmpty() ?
                            Collections.emptyListIterator() :
                            Collections.singletonList(doubleMSM(pairs)).iterator();
        }).reduce((e1, e2) -> new Tuple2<>(e1._1.add(e2._1), e1._2.add(e2._2)));
    }

    /* Used for profiling only */
    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT distributedSortedMSM(final JavaRDD<Tuple2<FieldT, GroupT>> input) {

        return input.mapPartitions(partition -> {
            final List<Tuple2<BigInteger, GroupT>> pairs = new ArrayList<>();

            while (partition.hasNext()) {
                final Tuple2<FieldT, GroupT> pair = partition.next();
                final BigInteger scalar = pair._1.toBigInteger();
                if (scalar.equals(BigInteger.ZERO)) {
                    continue;
                }
                pairs.add(new Tuple2<>(scalar, pair._2));
            }

            if (pairs.isEmpty()) {
                return Collections.emptyListIterator();
            }

            return Collections.singletonList(sortedMSM(pairs)).iterator();
        }).reduce(GroupT::add);
    }

    /* Used for profiling only */
    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT distributedBosCosterMSM(final JavaRDD<Tuple2<FieldT, GroupT>> input) {

        return input.mapPartitions(partition -> {
            final List<Tuple2<BigInteger, GroupT>> pairs = new ArrayList<>();

            while (partition.hasNext()) {
                Tuple2<FieldT, GroupT> part = partition.next();
                pairs.add(new Tuple2<>(part._1().toBigInteger(), part._2()));
            }

            if (pairs.isEmpty()) {
                return Collections.emptyListIterator();
            }

            return Collections.singletonList(bosCosterMSM(pairs)).iterator();
        }).reduce(GroupT::add);
    }

    /* Used for profiling only */
    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT distributedPippengerMSM(final JavaRDD<Tuple2<FieldT, GroupT>> input) {

        return input.mapPartitions(partition -> {
            final List<Tuple2<BigInteger, GroupT>> pairs = new ArrayList<>();

            int numBits = 0;
            while (partition.hasNext()) {
                final Tuple2<FieldT, GroupT> part = partition.next();
                final BigInteger scalar = part._1().toBigInteger();

                pairs.add(new Tuple2<>(scalar, part._2()));
                numBits = Math.max(numBits, scalar.bitLength());
            }

            if (pairs.isEmpty()) {
                return Collections.emptyListIterator();
            }

            return Collections.singletonList(pippengerMSM(pairs, numBits)).iterator();
        }).reduce(GroupT::add);
    }
}