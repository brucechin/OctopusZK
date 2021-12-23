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

    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    T serialMSMTest(
            final int scalarSize,
            final int windowSize,
            final List<List<T>> multiplesOfBase,
            final FieldT scalar) {

        final int outerc = (scalarSize + windowSize - 1) / windowSize;
        //System.out.println("JAVA out_len="+outerc + " inner_len="+ windowSize);
        final BigInteger bigScalar = scalar.toBigInteger();

        T res = multiplesOfBase.get(0).get(0);
        for (int outer = 0; outer < outerc; ++outer) {
            int inner = 0;
            for (int i = 0; i < windowSize; ++i) {
                if (bigScalar.testBit(outer * windowSize + i)) {         
                    inner |= 1 << i;
                }
            }
            //System.out.println("JAVA outer=" + outer + " inner=" +inner);
            res = res.add(multiplesOfBase.get(outer).get(inner));
            //System.out.println("target for add is=" + multiplesOfBase.get(outer).get(inner).toString());
        }
        //System.out.println("\n\n\n");

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
        
       if(multiplesOfBase.get(0).get(0).getClass().getName().equals("algebra.curves.barreto_naehrig.bn254a.BN254aG1")){
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
                    //System.out.println("out=" + i +" in=" + j + " X=" + byteToString(three_values.get(0).toByteArray()));
                    //System.out.println("three_values:" +three_values.get(0) + " " + three_values.get(1) + " " + three_values.get(2));
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
                bigScalars.add(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));    
            }
            long finish = System.currentTimeMillis();
            long timeElapsed = finish - start;
            System.out.println("data transfer preparation time elapsed: " + timeElapsed + " ms");
            byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, byteArrayX, byteArrayY, byteArrayZ, bigScalars, 1);

            start = System.currentTimeMillis();
            int size_of_bigint_cpp_side = 64;
            final List<T> jni_res = new ArrayList<>(scalars.size());
            //BigInteger G1_modulus = new BigInteger("21888242871839275222246405745257275088696311157297823662689037894645226208583");
            

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
                // bi_X = bi_X.mod(G1_modulus);
                // bi_Y = bi_Y.mod(G1_modulus);
                // bi_Z = bi_Z.mod(G1_modulus);

                T temp = multiplesOfBase.get(0).get(0).zero();
                temp.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
                jni_res.add(temp);
                //System.out.println("CUDA FixedBaseMSM output=" +temp.toString());
            }

            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
            return jni_res;
        }

        System.out.println("for BN254G2, we use the old way.");
        for (FieldT scalar : scalars) {
            T temp = serialMSM(scalarSize, windowSize, multiplesOfBase, scalar);
            res.add(temp);
            //System.out.println("JAVA FixedBaseMSM output=" + temp.toString());
        }
        //For BN254G2 we have not implemented it yet.
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
            final ArrayList<ArrayList<byte[]>> multiplesOfBase1_X,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase1_Y,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase1_Z,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Xa,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Ya,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Za,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Xb,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Yb,
            final ArrayList<ArrayList<byte[]>> multiplesOfBase2_Zb,
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

        ArrayList<ArrayList<byte[]>> multiplesOfBase1_X = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase1_Y = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase1_Z = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Xa = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Ya = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Za = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Xb = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Yb = new ArrayList<ArrayList<byte[]>>();
        ArrayList<ArrayList<byte[]>> multiplesOfBase2_Zb = new ArrayList<ArrayList<byte[]>>();

        int out_size1 = multiplesOfBase1.size();
        int in_size1 = multiplesOfBase1.get(0).size();
        int out_size2 = multiplesOfBase2.size();
        int in_size2 = multiplesOfBase2.get(0).size();
        for(int i =0; i < out_size1; i++){
            ArrayList<byte[]> tmpX = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpY = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpZ = new ArrayList<byte[]>();            
            for(int j = 0; j< in_size1; j++){
                ArrayList<BigInteger> three_values = multiplesOfBase1.get(i).get(j).BN254G1ToBigInteger();
                tmpX.add(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                tmpY.add(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                tmpZ.add(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
            }
            multiplesOfBase1_X.add(tmpX);
            multiplesOfBase1_Y.add(tmpY);
            multiplesOfBase1_Z.add(tmpZ);
        }
        for(int i =0; i < out_size2; i++){
            ArrayList<byte[]> tmpXa = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpYa = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpZa = new ArrayList<byte[]>();     
            ArrayList<byte[]> tmpXb = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpYb = new ArrayList<byte[]>();
            ArrayList<byte[]> tmpZb = new ArrayList<byte[]>();     
            for(int j = 0; j< in_size2; j++){
                ArrayList<BigInteger> six_values = multiplesOfBase2.get(i).get(j).BN254G2ToBigInteger();
                tmpXa.add(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                tmpXa.add(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                tmpYa.add(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                tmpYb.add(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                tmpZa.add(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                tmpZb.add(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
            }
            multiplesOfBase2_Xa.add(tmpXa);
            multiplesOfBase2_Xb.add(tmpXb);
            multiplesOfBase2_Ya.add(tmpYa);
            multiplesOfBase2_Yb.add(tmpYb);
            multiplesOfBase2_Za.add(tmpZa);
            multiplesOfBase2_Zb.add(tmpZb);

        }

        final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
        final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;

        ArrayList<byte[]> bigScalars = new ArrayList<byte[]>();
        for (FieldT scalar : scalars) {
            bigScalars.add(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));
        }

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

        // //because each G1 value takes up 3 BigIntegers, and each G2 takes up 6 BigIntegers.
        // for(int i = 0; i < scalars.size(); i++){
        //     byte[] slice1 = Arrays.copyOfRange(resultByteArray, 9*i*size_of_bigint_cpp_side, (9 * i + 3) * size_of_bigint_cpp_side);
        //     byte[] converted_back_X = new byte[64];
        //     byte[] converted_back_Y = new byte[64];
        //     byte[] converted_back_Z = new byte[64];
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_X[j] = slice1[size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Y[j] = slice1[2*size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Z[j] = slice1[3*size_of_bigint_cpp_side - j - 1];
        //     }

        //     BigInteger bi_X = new BigInteger(converted_back_X);
        //     BigInteger bi_Y = new BigInteger(converted_back_Y);
        //     BigInteger bi_Z = new BigInteger(converted_back_Z);

        //     G1T temp1 = multiplesOfBase1.get(0).get(0).zero();
        //     temp1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);

        //     byte[] slice2 = Arrays.copyOfRange(resultByteArray, (9*i +3)*size_of_bigint_cpp_side, (9*i+9)*size_of_bigint_cpp_side);
        //     byte[] converted_back2 = new byte[64];
 
        //     byte[] converted_back_Xa = new byte[64];
        //     byte[] converted_back_Ya = new byte[64];
        //     byte[] converted_back_Za = new byte[64];
        //     byte[] converted_back_Xb = new byte[64];
        //     byte[] converted_back_Yb = new byte[64];
        //     byte[] converted_back_Zb = new byte[64];
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Xa[j] = slice1[size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Xb[j] = slice1[2*size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Ya[j] = slice1[3*size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Yb[j] = slice1[4*size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Za[j] = slice1[5*size_of_bigint_cpp_side - j - 1];
        //     }
        //     for(int j =0; j < size_of_bigint_cpp_side; j++){
        //         converted_back_Zb[j] = slice1[6*size_of_bigint_cpp_side - j - 1];
        //     }


        //     BigInteger bi_Xa = new BigInteger(converted_back_Xa);
        //     BigInteger bi_Ya = new BigInteger(converted_back_Ya);
        //     BigInteger bi_Za = new BigInteger(converted_back_Za);
        //     BigInteger bi_Xb = new BigInteger(converted_back_Xb);
        //     BigInteger bi_Yb = new BigInteger(converted_back_Yb);
        //     BigInteger bi_Zb = new BigInteger(converted_back_Zb);

        //     G2T temp2 = multiplesOfBase2.get(0).get(0).zero();
        //     temp2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);
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

