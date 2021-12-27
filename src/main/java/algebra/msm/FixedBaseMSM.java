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

import org.apache.commons.collections.IteratorUtils;
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
import algebra.fields.Fp2;
import java.util.Iterator;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Arrays;
import java.io.*;
import org.apache.spark.SparkEnv;
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
        //System.out.println("AlgebraMSMFixedBaseMSM getWindowSize");
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
        System.out.println("numWindows="+ numWindows +"windowSize="+ windowSize );
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
            final int out_len, final int inner_len, final int batch_size, final int scalarSize,
            final byte[] multiplesOfBaseBN254XYZ,
            final byte[] bigScalarsArrays,
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
        //byte[] res = new byte[(temp.length + 3)/ 4 * 4];
        byte[] res = new byte[32]; //254-bit BN254 biginteger is enough to be saved in 32bytes.
        for(int i = 0; i < temp.length; i++){
            res[temp.length - i - 1] = temp[i];
        }
        
        return res;

    }


    //lianke: this batchMSM is for serial zkSNARK
    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<T> batchMSM(
            final int scalarSize,
            final int windowSize,
            final List<List<T>> multiplesOfBase,
            final List<FieldT> scalars) throws Exception{
        System.out.println("scalarSize len :" + scalarSize);
        System.out.println("windowSize len :" + windowSize);
        System.out.println("batchMSM len :" + scalars.size());
        System.out.println("multiplesOfBase len : " + multiplesOfBase.size() + " " + multiplesOfBase.get(0).size());
        System.out.println("multiplesOfBase type : " + multiplesOfBase.get(0).get(0).getClass().getName());
        System.out.println("scalars type : " + scalars.get(0).getClass().getName());
        
    
            int out_size = multiplesOfBase.size();
            int in_size = multiplesOfBase.get(0).size();
            long start = System.currentTimeMillis();
            ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );

            for(int i =0; i < out_size; i++){
                for(int j = 0; j< in_size; j++){
                    ArrayList<BigInteger> three_values = multiplesOfBase.get(i).get(j).BN254G1ToBigInteger();
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));

                }

            }
            byte[] baseByteArrayXYZ = outputStream.toByteArray();

            final int outerc = (scalarSize + windowSize - 1) / windowSize;
            
            ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );

            for (FieldT scalar : scalars) {
                bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));
            }
            byte[] bigScalarByteArray =  bigScalarStream.toByteArray();

            long finish = System.currentTimeMillis();
            long timeElapsed = finish - start;
            System.out.println("data transfer preparation time elapsed: " + timeElapsed + " ms");
            byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, scalars.size(), scalarSize, baseByteArrayXYZ, bigScalarByteArray, 1);
            
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
                //System.out.println("CUDA FixedBaseMSM output=" +temp.toString());
            }

            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
            return jni_res;
        

    }


    //this batchMSMPartition is for distributed zkSNARK, each executor could work on multiple partition and each partition map this method.
    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<Tuple2<Long, T>> batchMSMPartition(
            final int scalarSize,
            final int windowSize,
            final List<List<T>> multiplesOfBase,
            final List<Tuple2< Long, FieldT>> scalars) throws Exception {
        
                int out_size = multiplesOfBase.size();
            int in_size = multiplesOfBase.get(0).size();
            long start = System.currentTimeMillis();
            ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );
            //TODO compute the table on GPU side directly

            for(int i =0; i < out_size; i++){
                for(int j = 0; j< in_size; j++){
                    ArrayList<BigInteger> three_values = multiplesOfBase.get(i).get(j).BN254G1ToBigInteger();
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
                }

            }
            byte[] baseByteArrayXYZ = outputStream.toByteArray();

            final int outerc = (scalarSize + windowSize - 1) / windowSize;
            
            ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );

            for (Tuple2<Long, FieldT> scalar : scalars) {
                bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar._2.toBigInteger()));
            }
            byte[] bigScalarByteArray =  bigScalarStream.toByteArray();

            long finish = System.currentTimeMillis();
            long timeElapsed = finish - start;
            System.out.println("data transfer preparation time elapsed: " + timeElapsed + " ms");
            byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, scalars.size(),scalarSize,  baseByteArrayXYZ, bigScalarByteArray, 1);
            
            start = System.currentTimeMillis();
            int size_of_bigint_cpp_side = 64;
            final List<Tuple2<Long, T>> jni_res = new ArrayList<>(scalars.size());
            

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
                jni_res.add(new Tuple2<>(scalars.get(i)._1, temp));
                //System.out.println("CUDA FixedBaseMSM output=" +temp.toString());
            }

            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
            return jni_res;
        
    }


    public static <GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, GroupT>
    distributedBatchMSM(
            final int scalarSize,
            final int windowSize,
            final List<List<GroupT>> multiplesOfBase,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {
        long start1 = System.currentTimeMillis();

        final Broadcast<List<List<GroupT>>> baseBroadcast = sc.broadcast(multiplesOfBase);
        long finish1 = System.currentTimeMillis();
        long timeElapsed1 = finish1 - start1;
        System.out.println("Spark broadcast  time elapsed: " + timeElapsed1 + " ms");

        long start = System.currentTimeMillis();
        JavaPairRDD<Long, GroupT> true_result = scalars.mapPartitionsToPair(
            partition ->  {
                List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                List<Tuple2<Long, GroupT>> res =  batchMSMPartition(scalarSize, windowSize, baseBroadcast.getValue(), scalar_partition);
                return res.iterator();
            }
        );
        true_result.count();
        long finish = System.currentTimeMillis();
        long timeElapsed = finish - start;
        System.out.println("Spark processing time elapsed: " + timeElapsed + " ms");
        return true_result;
        // return scalars.mapToPair(scalar -> new Tuple2<>(
        //         scalar._1,
        //         serialMSM(scalarSize, windowSize, baseBroadcast.getValue(), scalar._2)));
    }
    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] doubleBatchMSMNativeHelper(
            final int outerc1,
            final int windowSize1,
            final int outerc2,
            final int windowSize2,
            final int out_len1, final int inner_len1, 
            final int out_len2, final int inner_len2, 
            final int batch_size,
            final byte[] multiplesOfBaseBN254G1XYZ,
            final byte[] multiplesOfBaseBN254G2XYZABC,
            final byte[] bigScalarsArrays);

    public static <G1T extends AbstractGroup<G1T>,
            G2T extends AbstractGroup<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> List<Tuple2<G1T, G2T>> doubleBatchMSM(
            final int scalarSize1,
            final int windowSize1,
            final List<List<G1T>> multiplesOfBase1,
            final int scalarSize2,
            final int windowSize2,
            final List<List<G2T>> multiplesOfBase2,
            final List<FieldT> scalars) throws Exception {

        final List<Tuple2<G1T, G2T>> res = new ArrayList<>(scalars.size());

        System.out.println("doubleBatchMSM info:");
        System.out.println("batchMSM size :" + scalars.size());
        System.out.println("scalarSize len :" + scalarSize1 + " " + scalarSize2 );
        System.out.println("windowSize len :" + windowSize1 + " " + windowSize2);
        System.out.println("multiplesOfBase1 len : " + multiplesOfBase1.size() + " " + multiplesOfBase1.get(0).size());
        System.out.println("multiplesOfBase2 len : " + multiplesOfBase2.size() + " " + multiplesOfBase2.get(0).size());
        System.out.println("multiplesOfBase1 type : " + multiplesOfBase1.get(0).get(0).getClass().getName());
        System.out.println("multiplesOfBase2 type : " + multiplesOfBase2.get(0).get(0).getClass().getName());
        System.out.println("scalars type : " + scalars.get(0).getClass().getName());

       
        int out_size1 = multiplesOfBase1.size();
        int in_size1 = multiplesOfBase1.get(0).size();
        int out_size2 = multiplesOfBase2.size();
        int in_size2 = multiplesOfBase2.get(0).size();
        long start = System.currentTimeMillis();

        ByteArrayOutputStream G1outputStream = new ByteArrayOutputStream( );

        ByteArrayOutputStream G2outputStream = new ByteArrayOutputStream( );
            //TODO compute the table on GPU side directly


        for(int i =0; i < out_size1; i++){
        
            for(int j = 0; j< in_size1; j++){
                ArrayList<BigInteger> three_values = multiplesOfBase1.get(i).get(j).BN254G1ToBigInteger();
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
            }

        }
        byte[] baseByteArrayXYZ = G1outputStream.toByteArray();

        for(int i =0; i < out_size2; i++){
   
            for(int j = 0; j< in_size2; j++){
                ArrayList<BigInteger> six_values = multiplesOfBase2.get(i).get(j).BN254G2ToBigInteger();
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
            }

        }
        byte[] baseByteArrayXYZABC = G2outputStream.toByteArray();

        final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
        final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;
        ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        for (FieldT scalar : scalars) {
            bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar.toBigInteger()));
        }
        long finish = System.currentTimeMillis();
        byte[] bigScalarByteArray =  bigScalarStream.toByteArray();

        long timeElapsed = finish - start;
        System.out.println("JAVA side prepare data time elapsed: " + timeElapsed + " ms");

        byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2,
                                        out_size1, in_size1, out_size2, in_size2, scalars.size(),
                                        baseByteArrayXYZ, baseByteArrayXYZABC, bigScalarByteArray);


        final List<Tuple2<G1T, G2T>> jni_res = new ArrayList<>(scalars.size());

        start = System.currentTimeMillis();
        int size_of_bigint_cpp_side = 64;

        // because each G1 value takes up 3 BigIntegers, and each G2 takes up 6 BigIntegers.
        for(int i = 0; i < scalars.size(); i++){
            byte[] slice1 = Arrays.copyOfRange(resultByteArray, 9*i*size_of_bigint_cpp_side, (9 * i + 3) * size_of_bigint_cpp_side);
            byte[] converted_back_X = new byte[64];
            byte[] converted_back_Y = new byte[64];
            byte[] converted_back_Z = new byte[64];
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_X[j] = slice1[size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Y[j] = slice1[2*size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Z[j] = slice1[3*size_of_bigint_cpp_side - j - 1];
            }

            BigInteger bi_X = new BigInteger(converted_back_X);
            BigInteger bi_Y = new BigInteger(converted_back_Y);
            BigInteger bi_Z = new BigInteger(converted_back_Z);
            //System.out.println("G1 X,Y,Z=" + bi_X + " " + bi_Y + " " + bi_Z);
            G1T temp1 = multiplesOfBase1.get(0).get(0).zero();
            temp1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);

            byte[] slice2 = Arrays.copyOfRange(resultByteArray, (9*i +3)*size_of_bigint_cpp_side, (9*i+9)*size_of_bigint_cpp_side);
 
            byte[] converted_back_Xa = new byte[64];
            byte[] converted_back_Ya = new byte[64];
            byte[] converted_back_Za = new byte[64];
            byte[] converted_back_Xb = new byte[64];
            byte[] converted_back_Yb = new byte[64];
            byte[] converted_back_Zb = new byte[64];
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Xa[j] = slice2[size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Xb[j] = slice2[2*size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Ya[j] = slice2[3*size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Yb[j] = slice2[4*size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Za[j] = slice2[5*size_of_bigint_cpp_side - j - 1];
            }
            for(int j =0; j < size_of_bigint_cpp_side; j++){
                converted_back_Zb[j] = slice2[6*size_of_bigint_cpp_side - j - 1];
            }


            BigInteger bi_Xa = new BigInteger(converted_back_Xa);
            BigInteger bi_Ya = new BigInteger(converted_back_Ya);
            BigInteger bi_Za = new BigInteger(converted_back_Za);
            BigInteger bi_Xb = new BigInteger(converted_back_Xb);
            BigInteger bi_Yb = new BigInteger(converted_back_Yb);
            BigInteger bi_Zb = new BigInteger(converted_back_Zb);

            G2T temp2 = multiplesOfBase2.get(0).get(0).zero();
            temp2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);
            //System.out.println("CUDA G2=" +temp2.toString());

            jni_res.add(new Tuple2<>(temp1, temp2));
        }

        finish = System.currentTimeMillis();
        timeElapsed = finish - start;
        System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");



        return jni_res;
    }

    public static <G1T extends AbstractGroup<G1T>,
        G2T extends AbstractGroup<G2T>,
        FieldT extends AbstractFieldElementExpanded<FieldT>> List<Tuple2< Long, Tuple2<G1T, G2T>>> doubleBatchMSMPartition(
        final int scalarSize1,
        final int windowSize1,
        final List<List<G1T>> multiplesOfBase1,
        final int scalarSize2,
        final int windowSize2,
        final List<List<G2T>> multiplesOfBase2,
        final List<Tuple2<Long, FieldT>> scalars) throws Exception{


            int out_size1 = multiplesOfBase1.size();
            int in_size1 = multiplesOfBase1.get(0).size();
            int out_size2 = multiplesOfBase2.size();
            int in_size2 = multiplesOfBase2.get(0).size();
            long start = System.currentTimeMillis();
    
            ByteArrayOutputStream G1outputStream = new ByteArrayOutputStream( );
    
            ByteArrayOutputStream G2outputStream = new ByteArrayOutputStream( );
    
            //TODO compute the table on GPU side directly
            for(int i =0; i < out_size1; i++){
            
                for(int j = 0; j< in_size1; j++){
                    ArrayList<BigInteger> three_values = multiplesOfBase1.get(i).get(j).BN254G1ToBigInteger();
                    G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                    G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                    G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
                }
    
            }
            byte[] baseByteArrayXYZ = G1outputStream.toByteArray();
    
            for(int i =0; i < out_size2; i++){
       
                for(int j = 0; j< in_size2; j++){
                    ArrayList<BigInteger> six_values = multiplesOfBase2.get(i).get(j).BN254G2ToBigInteger();
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                    G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
                }
    
            }
            byte[] baseByteArrayXYZABC = G2outputStream.toByteArray();
    
            final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
            final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;
            ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
            for (Tuple2<Long, FieldT> scalar : scalars) {
                bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar._2.toBigInteger()));
            }
            long finish = System.currentTimeMillis();
            byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
    
            long timeElapsed = finish - start;
            System.out.println("JAVA side prepare data time elapsed: " + timeElapsed + " ms");
    
            byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2,
                                            out_size1, in_size1, out_size2, in_size2, scalars.size(),
                                            baseByteArrayXYZ, baseByteArrayXYZABC, bigScalarByteArray);
    
    
            final List<Tuple2< Long, Tuple2<G1T, G2T>>> jni_res = new ArrayList<>(scalars.size());
    
            start = System.currentTimeMillis();
            int size_of_bigint_cpp_side = 64;
    
            // because each G1 value takes up 3 BigIntegers, and each G2 takes up 6 BigIntegers.
            for(int i = 0; i < scalars.size(); i++){
                byte[] slice1 = Arrays.copyOfRange(resultByteArray, 9*i*size_of_bigint_cpp_side, (9 * i + 3) * size_of_bigint_cpp_side);
                byte[] converted_back_X = new byte[64];
                byte[] converted_back_Y = new byte[64];
                byte[] converted_back_Z = new byte[64];
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_X[j] = slice1[size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Y[j] = slice1[2*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Z[j] = slice1[3*size_of_bigint_cpp_side - j - 1];
                }
    
                BigInteger bi_X = new BigInteger(converted_back_X);
                BigInteger bi_Y = new BigInteger(converted_back_Y);
                BigInteger bi_Z = new BigInteger(converted_back_Z);
                //System.out.println("G1 X,Y,Z=" + bi_X + " " + bi_Y + " " + bi_Z);
                G1T temp1 = multiplesOfBase1.get(0).get(0).zero();
                temp1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
    
                byte[] slice2 = Arrays.copyOfRange(resultByteArray, (9*i +3)*size_of_bigint_cpp_side, (9*i+9)*size_of_bigint_cpp_side);
     
                byte[] converted_back_Xa = new byte[64];
                byte[] converted_back_Ya = new byte[64];
                byte[] converted_back_Za = new byte[64];
                byte[] converted_back_Xb = new byte[64];
                byte[] converted_back_Yb = new byte[64];
                byte[] converted_back_Zb = new byte[64];
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Xa[j] = slice2[size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Xb[j] = slice2[2*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Ya[j] = slice2[3*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Yb[j] = slice2[4*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Za[j] = slice2[5*size_of_bigint_cpp_side - j - 1];
                }
                for(int j =0; j < size_of_bigint_cpp_side; j++){
                    converted_back_Zb[j] = slice2[6*size_of_bigint_cpp_side - j - 1];
                }
    
    
                BigInteger bi_Xa = new BigInteger(converted_back_Xa);
                BigInteger bi_Ya = new BigInteger(converted_back_Ya);
                BigInteger bi_Za = new BigInteger(converted_back_Za);
                BigInteger bi_Xb = new BigInteger(converted_back_Xb);
                BigInteger bi_Yb = new BigInteger(converted_back_Yb);
                BigInteger bi_Zb = new BigInteger(converted_back_Zb);
    
                G2T temp2 = multiplesOfBase2.get(0).get(0).zero();
                temp2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);
                //System.out.println("CUDA G2=" +temp2.toString());
    
                jni_res.add(new Tuple2<>(scalars.get(i)._1, new Tuple2<>(temp1, temp2)));
            }
    
            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
    
    
    
            return jni_res;

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
        return scalars.mapPartitionsToPair(
            partition ->  {
                    List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                    List<Tuple2<Long, Tuple2<G1T, G2T>>> res =  doubleBatchMSMPartition(scalarSize1, windowSize1, baseBroadcast1.getValue(), 
                                                                scalarSize2, windowSize2, baseBroadcast2.getValue(), scalar_partition);
                    return res.iterator();
            }
        );
        // return scalars.mapToPair(scalar -> new Tuple2<>(
        //         scalar._1,
        //         new Tuple2<>(
        //                 serialMSM(scalarSize1, windowSize1, baseBroadcast1.value(), scalar._2),
        //                 serialMSM(scalarSize2, windowSize2, baseBroadcast2.value(), scalar._2))));
    }
}

