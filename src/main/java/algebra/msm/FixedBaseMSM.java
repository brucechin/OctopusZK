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
            final ArrayList<ArrayList<byte[]>> multiplesOfBaseBN254X,
            final ArrayList<ArrayList<byte[]>> multiplesOfBaseBN254Y,
            final ArrayList<ArrayList<byte[]>> multiplesOfBaseBN254Z,
            final ArrayList<byte[]> bigScalars,
            int BNType);//1 is G1, 2 is G2.


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

    public static byte[] bigIntegerToByteArrayHelperCGBN(BigInteger bigint){
        byte[] temp = bigint.toByteArray();

        byte[] res = new byte[(temp.length + 3)/ 4 * 4];
        int new_len = (temp.length + 3) / 4 * 4;
        for(int i = 0; i < temp.length; i++){
            res[temp.length - i - 1] = temp[i];
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
        System.out.println("scalarSize len :" + scalarSize);
        System.out.println("windowSize len :" + windowSize);

        System.out.println("batchMSM len :" + scalars.size());
        System.out.println("multiplesOfBase len : " + multiplesOfBase.size() + " " + multiplesOfBase.get(0).size());
        System.out.println("multiplesOfBase type : " + multiplesOfBase.get(0).get(0).getClass().getName());
        System.out.println("scalars type : " + scalars.get(0).getClass().getName());
        
       // if(multiplesOfBase.get(0).get(0).getClass().getName().equals("algebra.curves.barreto_naehrig.bn254a.BN254aG1")){
        if(false){
        ArrayList<ArrayList<byte[]>> byteArrayX = new ArrayList<ArrayList<byte[]>>();
            ArrayList<ArrayList<byte[]>> byteArrayY = new ArrayList<ArrayList<byte[]>>();
            ArrayList<ArrayList<byte[]>> byteArrayZ = new ArrayList<ArrayList<byte[]>>();
    
            int out_size = multiplesOfBase.size();
            int in_size = multiplesOfBase.get(0).size();
            long start = System.currentTimeMillis();
    
            for(int i =0; i < out_size; i++){
                ArrayList<byte[]> tmpX = new ArrayList<byte[]>();
                ArrayList<byte[]> tmpY = new ArrayList<byte[]>();
                ArrayList<byte[]> tmpZ = new ArrayList<byte[]>();
                for(int j = 0; j< in_size; j++){
                    ArrayList<BigInteger> three_values = multiplesOfBase.get(i).get(j).BN254G1ToBigInteger();
                    tmpX.add(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                    tmpY.add(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                    tmpZ.add(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
                }
                byteArrayX.add(tmpX);
                byteArrayY.add(tmpY);
                byteArrayZ.add(tmpZ);
            }
    
            final int outerc = (scalarSize + windowSize - 1) / windowSize;
            ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
            for (FieldT scalar : scalars) {
                bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));    
            }
            long finish = System.currentTimeMillis();
            long timeElapsed = finish - start;
            System.out.println("data transfer preparation time elapsed: " + timeElapsed + " ms");
            byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, byteArrayX, byteArrayY, byteArrayZ, bigScalars, 1);

            start = System.currentTimeMillis();
            int size_of_bigint_cpp_side = 64;
            final List<T> jni_res = new ArrayList<>(scalars.size());

            for(int i = 0; i < scalars.size(); i++){
                byte[] slice = Arrays.copyOfRange(resultByteArray, 3*i*size_of_bigint_cpp_side, 3*(i+1)*size_of_bigint_cpp_side);

                byte[] converted_back_X = new byte[64];
                byte[] converted_back_Y = new byte[64];
                byte[] converted_back_Z = new byte[64];

                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_X[j] = slice[size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Y[j] = slice[2*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Z[j] = slice[3*size_of_bigint_cpp_side - j - 1];
                }

                BigInteger bi_X = new BigInteger(converted_back_X);
                BigInteger bi_Y = new BigInteger(converted_back_Y);
                BigInteger bi_Z = new BigInteger(converted_back_Z);

                T temp = multiplesOfBase.get(0).get(0).zero();
                temp.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
                jni_res.add(temp);
            
            }

            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
            return jni_res;
        }else{
            System.out.println("for BN254G2, we use the old way.");
            for (FieldT scalar : scalars) {
                res.add(serialMSM(scalarSize, windowSize, multiplesOfBase, scalar));
            }
            //For BN254G2 we have not implemented it yet.
            return res;
        }
        


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
        System.out.println("scalars type : " + scalars.get(0).getClass().getName());

        ArrayList<ArrayList<byte[]>> byteArray1 = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> byteArray2 = new ArrayList<ArrayList<byte[]>>();

        // int out_size1 = multiplesOfBase1.size();
        // int in_size1 = multiplesOfBase1.get(0).size();
        // int out_size2 = multiplesOfBase2.size();
        // int in_size2 = multiplesOfBase2.get(0).size();
        // //TODO lianke : parallelize it
        // for(int i =0; i < out_size1; i++){
        //     ArrayList<byte[]> tmp = new ArrayList<byte[]>();
        //     for(int j = 0; j< in_size1; j++){
        //         tmp.add(bigIntegerToByteArrayHelper(multiplesOfBase1.get(i).get(j).toBigInteger()));
        //     }
        //     byteArray1.add(tmp);
        // }
        // for(int i =0; i < out_size2; i++){
        //     ArrayList<byte[]> tmp = new ArrayList<byte[]>();
        //     for(int j = 0; j< in_size2; j++){
        //         tmp.add(bigIntegerToByteArrayHelper(multiplesOfBase2.get(i).get(j).toBigInteger()));
        //     }
        //     byteArray2.add(tmp);
        // }

        // final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
        // final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;

        // ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        // for (FieldT scalar : scalars) {
        //     bigScalars.add(bigIntegerToByteArrayHelper(scalar.toBigInteger()));
        // }

        //TODO lianke for verification purpose
        for (FieldT scalar : scalars) {
            res.add(new Tuple2<>(
                    serialMSM(scalarSize1, windowSize1, multiplesOfBase1, scalar),
                    serialMSM(scalarSize2, windowSize2, multiplesOfBase2, scalar)));
        }


        // byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2, byteArray1, byteArray2, bigScalars);


        // final List<Tuple2<G1T, G2T>> jni_res = new ArrayList<>(scalars.size());

        // long start = System.currentTimeMillis();
        // int size_of_bigint_cpp_side = 64;
        // BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        
        // for(int i = 0; i < scalars.size(); i++){
        //     byte[] slice1 = Arrays.copyOfRange(resultByteArray, 2*i*size_of_bigint_cpp_side, (2*i+1)*size_of_bigint_cpp_side);
        //     byte[] converted_back1 = new byte[64];
        //     for(int j = 63; j >= 3; j-=4){
        //         converted_back1[j] = slice1[j - 3];
        //         converted_back1[j-1] = slice1[j - 2];
        //         converted_back1[j-2] = slice1[j - 1];
        //         converted_back1[j-3] = slice1[j ];
        //     }
        //     BigInteger bi1 = new BigInteger(converted_back1);
        //     BigInteger output1 = bi1.mod(modulus);
        //     G1T temp1 = multiplesOfBase1.get(0).get(0).zero();
        //     temp1.setBigInteger(output1);

        //     byte[] slice2 = Arrays.copyOfRange(resultByteArray, (2*i +1)*size_of_bigint_cpp_side, (2*i+2)*size_of_bigint_cpp_side);
        //     byte[] converted_back2 = new byte[64];
        //     for(int j = 63; j >= 3; j-=4){
        //         converted_back2[j] = slice2[j - 3];
        //         converted_back2[j-1] = slice2[j - 2];
        //         converted_back2[j-2] = slice2[j - 1];
        //         converted_back2[j-3] = slice2[j ];
        //     }
        //     BigInteger bi2 = new BigInteger(converted_back2);
        //     BigInteger output2 = bi2.mod(modulus);
        //     G2T temp2 = multiplesOfBase2.get(0).get(0).zero();
        //     temp2.setBigInteger(output2);
        //     jni_res.add(new Tuple2<>(temp1, temp2));
        // }

        // long finish = System.currentTimeMillis();
        // long timeElapsed = finish - start;
        // System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");




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

