/* @file
 *****************************************************************************
 * @author     This file is part of zkspark, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

package algebra.msm;
import java.lang.reflect.*;
import algebra.fields.AbstractFieldElementExpanded;
import algebra.groups.AbstractGroup;
import algebra.curves.AbstractG1;
import algebra.curves.fake.FakeG1;
import algebra.curves.fake.FakeG2;
import algebra.curves.fake.fake_parameters.FakeFqParameters;
import algebra.curves.AbstractG2;
import org.apache.spark.api.java.JavaPairRDD;
import org.apache.spark.api.java.JavaSparkContext;
//import com.ibm.gpuenabler.*;
import org.apache.spark.broadcast.Broadcast;
import scala.Tuple2;
import algebra.fields.Fp;
import algebra.fields.fieldparameters.LargeFpParameters;
import algebra.groups.AdditiveIntegerGroup;
import algebra.groups.integergroupparameters.LargeAdditiveIntegerGroupParameters;
import java.math.BigInteger;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Arrays;
import java.io.*;

public class FixedBaseMSM {

    /**
     * Compute table of window sizes.
     */
    static {
		System.loadLibrary("AlgebraMSMFixedBaseMSM");
        System.out.println("AlgebraMSMFixedBaseMSM loaded");
	}

    // public static void main(String[] args) {

    //     System.out.println("FixedBaseMSM java side started");
    //     LargeAdditiveIntegerGroupParameters GroupParameters = new LargeAdditiveIntegerGroupParameters();
    //     final AdditiveIntegerGroup base = new AdditiveIntegerGroup(7, GroupParameters);
    //     final LargeFpParameters FpParameters = new LargeFpParameters();
    //     final Fp scalar = new Fp("20000000000", FpParameters);

    //     final int scalarSize = scalar.bitSize();
    //     final int windowSize = 2;

    //     List<List<AdditiveIntegerGroup>> windowTable = FixedBaseMSM
    //             .getWindowTable(base, scalarSize, windowSize);
    //     AdditiveIntegerGroup result = FixedBaseMSM
    //             .serialMSM(scalarSize, windowSize, windowTable, scalar);
    //     AdditiveIntegerGroup answers = new AdditiveIntegerGroup(1400, GroupParameters);
    //     System.out.println("FixedBaseMSM java side finished");

    // }
    public static <GroupT extends AbstractGroup<GroupT>> int getWindowSize(
            final long numScalars, final GroupT groupFactory) {

        if (groupFactory.fixedBaseWindowTable().isEmpty()) {
            return 17;
        }

        long window = 1;
        for (int i = groupFactory.fixedBaseWindowTable().size() - 1; i >= 0; i--) {
            final int value = groupFactory.fixedBaseWindowTable().get(i);
            if (value != 0 && numScalars >= value) {
                window = i + 1;
                break;
            }
        }

        return window > Integer.MAX_VALUE ? Integer.MAX_VALUE : (int) window;
    }

    /**
     * Computes the window table for a given base element.
     */
    public static <GroupT extends AbstractGroup<GroupT>> List<List<GroupT>> getWindowTable(
            final GroupT base, final int scalarSize, final int windowSize) {
        final int numWindows =
                (scalarSize % windowSize == 0) ? scalarSize / windowSize : scalarSize / windowSize + 1;
        final int innerLimit = (int) Math.pow(2, windowSize);

        final List<List<GroupT>> windowTable = new ArrayList<>();

        // If window table size is 0, return just the zero element.
        if (numWindows == 0) {
            return Collections.singletonList(Collections.singletonList(base.zero()));
        }

        GroupT baseOuter = base;
        for (int outer = 0; outer < numWindows; outer++) {
            windowTable.add(new ArrayList<>(innerLimit));
            GroupT baseInner = base.zero();
            for (int inner = 0; inner < innerLimit; inner++) {
                windowTable.get(outer).add(baseInner);
                baseInner = baseInner.add(baseOuter);
            }

            for (int w = 0; w < windowSize; w++) {
                baseOuter = baseOuter.twice();
            }
        }

        return windowTable;
    }

    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] batchMSMNativeHelper(
            final int outerc,
            final int windowSize,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase,
            final ArrayList<byte[]> bigScalars);

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




    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    T serialMSM(
            final int scalarSize,
            final int windowSize,
            final List<List<T>> multiplesOfBase,
            final FieldT scalar) {

        final int outerc = (scalarSize + windowSize - 1) / windowSize;
        final BigInteger bigScalar = scalar.toBigInteger();

        T res = multiplesOfBase.get(0).get(0);

        for (int outer = 0; outer < outerc; ++outer) {
            int inner = 0;
            for (int i = 0; i < windowSize; ++i) {
                if (bigScalar.testBit(outer * windowSize + i)) {
                    
                    inner |= 1 << i;
                }
            }

            res = res.add(multiplesOfBase.get(outer).get(inner));

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

    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<T> batchMSM(
            final int scalarSize,
            final int windowSize,
            final List<List<T>> multiplesOfBase,
            final List<FieldT> scalars) {
        final List<T> res = new ArrayList<>(scalars.size());
        //System.out.println("scalarSize len :" + scalarSize);
        //System.out.println("windowSize len :" + windowSize);

        System.out.println("batchMSM len :" + scalars.size());
        //TODO lianke here if batchMSM len is smaller than some threshold, we can do it on CPU. Here, batchMSM and multiplesOfBase size should be large(somehow equal to the number of constraints).
        System.out.println("multiplesOfBase len : " + multiplesOfBase.size() + " " + multiplesOfBase.get(0).size());
        System.out.println("multiplesOfBase type : " + multiplesOfBase.get(0).get(0).getClass().getName());
        
        ArrayList<ArrayList<byte[]>> byteArray = new ArrayList<ArrayList<byte[]>>();
        int out_size = multiplesOfBase.size();
        int in_size = multiplesOfBase.get(0).size();
        long start = System.currentTimeMillis();

        for(int i =0; i < out_size; i++){
            ArrayList<byte[]> tmp = new ArrayList<byte[]>();
            for(int j = 0; j< in_size; j++){
                tmp.add(bigIntegerToByteArrayHelper(multiplesOfBase.get(i).get(j).toBigInteger()));
                //System.out.println(multiplesOfBase.get(i).get(j).toBigInteger().toByteArray().length);
            }
            byteArray.add(tmp);
        }

        final int outerc = (scalarSize + windowSize - 1) / windowSize;
        ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        for (FieldT scalar : scalars) {

            bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));
            if(scalars.size() == 4){
                System.out.println(byteToString(scalar.toBigInteger().toByteArray()));
            }
        }
        long finish = System.currentTimeMillis();
        long timeElapsed = finish - start;
        System.out.println("data transfer preparation time elapsed: " + timeElapsed + " ms");
        byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, byteArray, bigScalars);


        for (FieldT scalar : scalars) {
            res.add(serialMSM(scalarSize, windowSize, multiplesOfBase, scalar));
        }
        final List<T> jni_res = new ArrayList<>(scalars.size());
        //TODO move modulus to cpp side to accelerate more.



        

        start = System.currentTimeMillis();
        int size_of_bigint_cpp_side = 64;
        BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        for(int i = 0; i < scalars.size(); i++){
            byte[] slice = Arrays.copyOfRange(resultByteArray, i*size_of_bigint_cpp_side, (i+1)*size_of_bigint_cpp_side);//in cpp side, BigInt is 32 bytes.

            byte[] converted_back = new byte[64];
            for(int j = 63; j >= 3; j-=4){
                converted_back[j] = slice[j - 3];
                converted_back[j-1] = slice[j - 2];
                converted_back[j-2] = slice[j - 1];
                converted_back[j-3] = slice[j ];
            }

            BigInteger bi = new BigInteger(converted_back);
            BigInteger output = bi.mod(modulus);
            T temp = multiplesOfBase.get(0).get(0).zero();
            temp.setBigInteger(output);
            jni_res.add(temp);
            // if(i == 0){
            //     System.out.println("java side fixedbase msm slice " + byteToString(slice));
            //     System.out.println("java side fixedbase msm convertedback " + byteToString(converted_back));
            //     System.out.println("java side fixedbase msm bi" + byteToString(bi.toByteArray()));  
            //     System.out.println("java side fixedbase msm res " + byteToString(bigIntegerToByteArrayHelper(res.get(i).toBigInteger())));
            //     System.out.println("java side fixedbase msm jni res " + byteToString(bigIntegerToByteArrayHelper(temp.toBigInteger())));  
            // }
 

            if(!res.get(i).toBigInteger().equals(temp.toBigInteger())){
                System.out.println("error in FixedBaseMSM.batchMSM JNI computation");
            }
            // System.out.println(
            // "\n                 res:" + res.get(i).toBigInteger()+
            // "\n jni modulus output :" + temp.toBigInteger() + " " + res.get(i).toBigInteger().equals(temp.toBigInteger()) );
        
        }

        finish = System.currentTimeMillis();
        timeElapsed = finish - start;
        System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
        return res;
    }

    public static <GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, GroupT>
    distributedBatchMSM(
            final int scalarSize,
            final int windowSize,
            final List<List<GroupT>> multiplesOfBase,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {

        final Broadcast<List<List<GroupT>>> baseBroadcast = sc.broadcast(multiplesOfBase);

        return scalars.mapToPair(scalar -> new Tuple2<>(
                scalar._1,
                serialMSM(scalarSize, windowSize, baseBroadcast.getValue(), scalar._2)));
    }
    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] doubleBatchMSMNativeHelper(
            final int outerc1,
            final int windowSize1,
            final int outerc2,
            final int windowSize2,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase1,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2,
            final ArrayList<byte[]> bigScalars);

    public static <G1T extends AbstractGroup<G1T>,
            G2T extends AbstractGroup<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> List<Tuple2<G1T, G2T>> doubleBatchMSM(
            final int scalarSize1,
            final int windowSize1,
            final List<List<G1T>> multiplesOfBase1,
            final int scalarSize2,
            final int windowSize2,
            final List<List<G2T>> multiplesOfBase2,
            final List<FieldT> scalars) {

        final List<Tuple2<G1T, G2T>> res = new ArrayList<>(scalars.size());

        System.out.println("doubleBatchMSM info:");
        System.out.println("batchMSM size :" + scalars.size());
        System.out.println("scalarSize len :" + scalarSize1 + " " + scalarSize2 );
        System.out.println("windowSize len :" + windowSize1 + " " + windowSize2);
        //TODO lianke here if batchMSM len is smaller than some threshold, we can do it on CPU. Here, batchMSM and multiplesOfBase size should be large(somehow equal to the number of constraints).
        System.out.println("multiplesOfBase1 len : " + multiplesOfBase1.size() + " " + multiplesOfBase1.get(0).size());
        System.out.println("multiplesOfBase2 len : " + multiplesOfBase2.size() + " " + multiplesOfBase2.get(0).size());
        System.out.println("multiplesOfBase1 type : " + multiplesOfBase1.get(0).get(0).getClass().getName());
        System.out.println("multiplesOfBase2 type : " + multiplesOfBase2.get(0).get(0).getClass().getName());

        ArrayList<ArrayList<byte[]>> byteArray1 = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> byteArray2 = new ArrayList<ArrayList<byte[]>>();

        int out_size1 = multiplesOfBase1.size();
        int in_size1 = multiplesOfBase1.get(0).size();
        int out_size2 = multiplesOfBase2.size();
        int in_size2 = multiplesOfBase2.get(0).size();
        //TODO lianke : parallelize it
        for(int i =0; i < out_size1; i++){
            ArrayList<byte[]> tmp = new ArrayList<byte[]>();
            for(int j = 0; j< in_size1; j++){
                tmp.add(bigIntegerToByteArrayHelper(multiplesOfBase1.get(i).get(j).toBigInteger()));
            }
            byteArray1.add(tmp);
        }
        for(int i =0; i < out_size2; i++){
            ArrayList<byte[]> tmp = new ArrayList<byte[]>();
            for(int j = 0; j< in_size2; j++){
                tmp.add(bigIntegerToByteArrayHelper(multiplesOfBase2.get(i).get(j).toBigInteger()));
            }
            byteArray2.add(tmp);
        }

        final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
        final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;

        ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        for (FieldT scalar : scalars) {
            bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));
        }

        //TODO lianke for verification purpose
        for (FieldT scalar : scalars) {
            res.add(new Tuple2<>(
                    serialMSM(scalarSize1, windowSize1, multiplesOfBase1, scalar),
                    serialMSM(scalarSize2, windowSize2, multiplesOfBase2, scalar)));
        }


        byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2, byteArray1, byteArray2, bigScalars);
        //TODO collect the result from JNI. should be trivial, but brings some overhead.


        final List<Tuple2<G1T, G2T>> jni_res = new ArrayList<>(scalars.size());

        long start = System.currentTimeMillis();
        int size_of_bigint_cpp_side = 64;
        BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        
        for(int i = 0; i < scalars.size(); i++){
            byte[] slice1 = Arrays.copyOfRange(resultByteArray, 2*i*size_of_bigint_cpp_side, (2*i+1)*size_of_bigint_cpp_side);
            byte[] converted_back1 = new byte[64];
            for(int j = 63; j >= 3; j-=4){
                converted_back1[j] = slice1[j - 3];
                converted_back1[j-1] = slice1[j - 2];
                converted_back1[j-2] = slice1[j - 1];
                converted_back1[j-3] = slice1[j ];
            }
            BigInteger bi1 = new BigInteger(converted_back1);
            BigInteger output1 = bi1.mod(modulus);
            G1T temp1 = multiplesOfBase1.get(0).get(0).zero();
            temp1.setBigInteger(output1);

            byte[] slice2 = Arrays.copyOfRange(resultByteArray, (2*i +1)*size_of_bigint_cpp_side, (2*i+2)*size_of_bigint_cpp_side);
            byte[] converted_back2 = new byte[64];
            for(int j = 63; j >= 3; j-=4){
                converted_back2[j] = slice2[j - 3];
                converted_back2[j-1] = slice2[j - 2];
                converted_back2[j-2] = slice2[j - 1];
                converted_back2[j-3] = slice2[j ];
            }
            BigInteger bi2 = new BigInteger(converted_back2);
            BigInteger output2 = bi2.mod(modulus);
            G2T temp2 = multiplesOfBase2.get(0).get(0).zero();
            temp2.setBigInteger(output2);


            jni_res.add(new Tuple2<>(temp1, temp2));
            if(!res.get(i)._1.toBigInteger().equals(temp1.toBigInteger()) || !res.get(i)._2.toBigInteger().equals(temp2.toBigInteger())){
                System.out.println("error in FixedBaseMSM.doubleBatchMSM JNI computation");
            }
            // System.out.println(
            // "\n                 res:" + res.get(i).toBigInteger()+
            // "\n jni modulus output :" + temp.toBigInteger() + " " + res.get(i).toBigInteger().equals(temp.toBigInteger()) );
        
        }

        long finish = System.currentTimeMillis();
        long timeElapsed = finish - start;
        System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");




        return res;
    }

    public static <G1T extends AbstractG1<G1T>,
            G2T extends AbstractG2<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, Tuple2<G1T, G2T>>
    distributedDoubleBatchMSM(
            final int scalarSize1,
            final int windowSize1,
            final List<List<G1T>> multiplesOfBase1,
            final int scalarSize2,
            final int windowSize2,
            final List<List<G2T>> multiplesOfBase2,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {

        final Broadcast<List<List<G1T>>> baseBroadcast1 = sc.broadcast(multiplesOfBase1);
        final Broadcast<List<List<G2T>>> baseBroadcast2 = sc.broadcast(multiplesOfBase2);

        return scalars.mapToPair(scalar -> new Tuple2<>(
                scalar._1,
                new Tuple2<>(
                        serialMSM(scalarSize1, windowSize1, baseBroadcast1.value(), scalar._2),
                        serialMSM(scalarSize2, windowSize2, baseBroadcast2.value(), scalar._2))));
    }
}

// |00001110|10111110|11110011|10010100|11001101|01011010|11110111|10000011|00011000|11010100|01101011|01010011|00110011|01010110|00111110|11100001|11101011|10101001|00101111|01110001|01011111|00111000|11010101|
// 00000000011110110100010001101111|11011101111111111101110001011000|00010110011001000010001010100110|00101001000110010001011010100110|11111111010001001011100101011111|         00100110101000101110101011101101|

