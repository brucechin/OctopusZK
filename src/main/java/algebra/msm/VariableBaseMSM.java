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
        for (int k = numGroups - 1; k >= 0; k--) {
            if (k < numGroups - 1) {
                for (int i = 0; i < c; i++) {
                    result = result.twice();
                }
            }

            final ArrayList<GroupT> buckets = new ArrayList<>(bucketsModel);

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


        }

        return result;
    }

    
    public static <GroupT extends AbstractGroup<GroupT>> ArrayList<GroupT> pippengerMSMTestPurpose(
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
        System.out.println("java side length: " + length + "log2 length: " + log2Length + "c " + c + " numGroups " + numGroups + "numBuckets " +  numBuckets);
        ArrayList<GroupT> resultArray = new ArrayList<>();
        for (int k = numGroups - 1; k >= 0; k--) {
            if (k < numGroups - 1) {
                for (int i = 0; i < c; i++) {
                    result = result.twice();
                }
            }
            // if(k > numGroups - 20) {
            //     System.out.println("java side k="+k + "  after exp result is :" + byteToString(result.toBigInteger().toByteArray()));
            // }
            final ArrayList<GroupT> buckets = new ArrayList<>(bucketsModel);

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

                // if(k == 80 && i  < 10){
                //     System.out.println("java side i="+i + "  after testBit, id is :" + id);
                // }

                // Potentially use mixed addition here.
                buckets.set(id, buckets.get(id).add(input.get(i)._2));
            }
            // if(k == 80 ) {
            //     for(int i = 0; i < numBuckets; i++){
            //         System.out.println("java side buckets index " + i + ": " + byteToString( buckets.get(i).toBigInteger().toByteArray()));
            //     }
    
            // }

            GroupT runningSum = zero;

            for (int i = numBuckets - 1; i > 0; i--) {
                // Potentially use mixed addition here.
                // if(k == numGroups - 1 ) {
                //         System.out.println("java side k " + k + " buckets index " + i + ", runningSum is " + byteToString( runningSum.toBigInteger().toByteArray()));
                //         System.out.println("result is :" + byteToString( result.toBigInteger().toByteArray()));
                    
                // }
                runningSum = runningSum.add(buckets.get(i));
                result = result.add(runningSum);
            }
            resultArray.add(result);
            // if(k < 10) {
            // System.out.println("java side numGroups="+ numGroups + "    l="+k + "  after adding runningSum, result is :" + byteToString(result.toBigInteger().toByteArray()));

            // }

        }

        return resultArray;
    }

    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] variableBaseSerialMSMNativeHelper(
            final ArrayList<byte[]> bases,
            final ArrayList<byte[]> scalars);

    public static <
            GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    GroupT serialMSM(final List<FieldT> scalars, final List<GroupT> bases) {
        // System.out.println("variableBaseSerialMSM info:");
        // System.out.println("variableBaseSerialMSM base size :" + bases.size() + " type:" +bases.get(0).getClass().getName());
        // System.out.println("variableBaseSerialMSM scalars size :" + scalars.size() + " type:"+scalars.get(0).getClass().getName());

        assert (bases.size() == scalars.size());




        // ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        // for (FieldT scalar : scalars) {
        //     bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));
        // }
        // ArrayList<byte[]> basesArray = new ArrayList<byte[]>();
        // for (GroupT base : bases){
        //     basesArray.add(bigIntegerToByteArrayHelper(base.toBigInteger()));
        // }

        // byte[] resArray = variableBaseSerialMSMNativeHelper(basesArray, bigScalars);
        
        // int size_of_bigint_cpp_side = 64;
        // BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        
        // byte[] converted_back = new byte[size_of_bigint_cpp_side];
        // for(int j = 63; j >= 3; j-=4){
        //     converted_back[j] = resArray[ j - 3];
        //     converted_back[j-1] = resArray[j - 2];
        //     converted_back[j-2] = resArray[ j - 1];
        //     converted_back[j-3] = resArray[ j];
        // }
        // BigInteger bi = new BigInteger(converted_back);
        // BigInteger output = bi.mod(modulus);
        // GroupT res = bases.get(0).zero();
        // res.setBigInteger(output);




        // lianke: below is the original code used for verify the JNI cpp side code computation is correct

        final List<Tuple2<BigInteger, GroupT>> filteredInput = new ArrayList<>();

        GroupT acc = bases.get(0).zero();


        int numBits = 0;
        for (int i = 0; i < bases.size(); i++) {
            final BigInteger scalar = scalars.get(i).toBigInteger();
            if (scalar.equals(BigInteger.ZERO)) {
                continue;
            }

            final GroupT base = bases.get(i);

            if (scalar.equals(BigInteger.ONE)) {
                acc = acc.add(base);
            } else {
                filteredInput.add(new Tuple2<>(scalar, base));
                //System.out.println("java side filteredInput add <scalar,base> index:" + i + "\n"  +  byteToString(scalar.toByteArray()) + ",\n " + byteToString(base.toBigInteger().toByteArray()));

                numBits = Math.max(numBits, scalar.bitLength());
            }
        }

       System.out.println("java side filteredInput size : " + filteredInput.size() + " number of bits " + numBits + " current acc is " + byteToString(acc.toBigInteger().toByteArray()));

       if (!filteredInput.isEmpty()) {
            acc = acc.add(pippengerMSM(filteredInput, numBits));
        }
        //ArrayList<GroupT> resultArray = new ArrayList<>();
        // if (!filteredInput.isEmpty()) {
        //     //resultArray = pippengerMSMTestPurpose(filteredInput, numBits);
        //     acc = acc.add(pippengerMSM(filteredInput, numBits));
        //     //System.out.println("java side final acc :" + byteToString(acc.toBigInteger().toByteArray()));
        // }

        // if(! res.toBigInteger().equals(acc.toBigInteger())){
        //     System.out.println("found error in pippengerMSM JNI computation");
        // }
    


        return acc;
    }


    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] variableBaseDoubleMSMNativeHelper(
            final ArrayList<byte[]> bases1,
            final ArrayList<byte[]> bases2,
            final ArrayList<byte[]> scalars);

    public static <
            T1 extends AbstractGroup<T1>,
            T2 extends AbstractGroup<T2>,
            FieldT extends AbstractFieldElementExpanded<FieldT>>
    Tuple2<T1, T2> doubleMSM(final List<FieldT> scalars, final List<Tuple2<T1, T2>> bases) {
        // System.out.println("variableBaseDoubleMSM info:");
        // System.out.println("variableBaseDoubleMSM base size :" + bases.size() + " type:" +bases.get(0)._1.getClass().getName() + " " +bases.get(0)._2.getClass().getName()  );
        // System.out.println("variableBaseDoubleMSM scalars size :" + scalars.size() + " type:"+scalars.get(0).getClass().getName());

        assert (bases.size() == scalars.size());

        final int size = bases.size();
        assert (size > 0);



        // ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        // for (FieldT scalar : scalars) {
        //     bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));
        // }
        // ArrayList<byte[]> basesArray1 = new ArrayList<byte[]>();
        // ArrayList<byte[]> basesArray2 = new ArrayList<byte[]>();

        // for (Tuple2<T1, T2> base : bases){
        //     basesArray1.add(bigIntegerToByteArrayHelper(base._1.toBigInteger()));
        //     basesArray2.add(bigIntegerToByteArrayHelper(base._2.toBigInteger()));
        // }

        // byte[] resArray = variableBaseDoubleMSMNativeHelper(basesArray1, basesArray2, bigScalars);


        // int size_of_bigint_cpp_side = 64;
        // BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        // byte[] converted_back1 = new byte[size_of_bigint_cpp_side];
        // for(int j = 63; j >= 3; j-=4){
        //     converted_back1[j] = resArray[j - 3];
        //     converted_back1[j-1] = resArray[j - 2];
        //     converted_back1[j-2] = resArray[j - 1];
        //     converted_back1[j-3] = resArray[j];
        // }
        // BigInteger bi = new BigInteger(converted_back1);
        // BigInteger output1 = bi.mod(modulus);
        // T1 res1 = bases.get(0)._1.zero();
        // res1.setBigInteger(output1);

        // byte[] converted_back2 = new byte[size_of_bigint_cpp_side];
        // for(int j = 63; j >= 3; j-=4){
        //     converted_back2[j] = resArray[64 + j - 3];
        //     converted_back2[j-1] = resArray[64 + j - 2];
        //     converted_back2[j-2] = resArray[64 + j - 1];
        //     converted_back2[j-3] = resArray[64 + j];
        // }
        // BigInteger bi2 = new BigInteger(converted_back2);
        // BigInteger output2 = bi2.mod(modulus);
        // T2 res2 = bases.get(0)._2.zero();
        // res2.setBigInteger(output2);





        // lianke: below is the original code used for verify the JNI cpp side code computation is correct

        final ArrayList<Tuple2<BigInteger, T1>> converted1 = new ArrayList<>(size);
        final ArrayList<Tuple2<BigInteger, T2>> converted2 = new ArrayList<>(size);

        T1 acc1 = bases.get(0)._1.zero();
        T2 acc2 = bases.get(0)._2.zero();
        int numBits = 0;

        for (int i = 0; i < size; i++) {
            final Tuple2<T1, T2> value = bases.get(i);
            final BigInteger scalar = scalars.get(i).toBigInteger();
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


        T1 true_res1 = converted1.isEmpty() ? acc1 : acc1.add(VariableBaseMSM.pippengerMSM(converted1, numBits));
        T2 true_res2 =  converted2.isEmpty() ? acc2 : acc2.add(VariableBaseMSM.pippengerMSM(converted2, numBits));
       
        // System.out.println("java side final acc1 :" + byteToString(true_res1.toBigInteger().toByteArray()));
        // System.out.println("java side final acc2 :" + byteToString(true_res2.toBigInteger().toByteArray()));

        // if(!true_res1.toBigInteger().equals(res1.toBigInteger())){
        //     System.out.println("error in first VariableBaseMSM .doubleMSM JNI computation");
        //     System.out.println(true_res1.toBigInteger() + " " + res1.toBigInteger());
        // }

        // if(!true_res2.toBigInteger().equals(res2.toBigInteger())){
        //     System.out.println("error in second VariableBaseMSM .doubleMSM JNI computation");
        //     System.out.println(true_res2.toBigInteger() + " " + res2.toBigInteger());
        // }
        return new Tuple2<>(true_res1, true_res2);




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