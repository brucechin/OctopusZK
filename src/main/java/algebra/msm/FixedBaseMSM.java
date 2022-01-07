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
import org.apache.spark.TaskContext;

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
            int BNType, int taskID);//1 is G1, 2 is G2.

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



    public static byte[] bigIntegerToByteArrayHelperCGBN(BigInteger bigint){
        byte[] temp = bigint.toByteArray();
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
            final int out_size, final int in_size,
            final T baseElement,
            final List<FieldT> scalars) throws Exception{
        // System.out.println("scalarSize len :" + scalarSize);
        // System.out.println("windowSize len :" + windowSize);
        // System.out.println("batchMSM len :" + scalars.size());
        // System.out.println("scalars type : " + scalars.get(0).getClass().getName());
        // System.out.println("base type :" + baseElement.getClass().getName());
    
        if(baseElement.getClass().getName() == "algebra.curves.barreto_naehrig.bn254a.BN254aG1"){
            //G1
            int G1_iteration_batch_size = (1 << 23);
            final List<T> jni_res = new ArrayList<>(scalars.size());
            //because when size is too large, we will exceed the capacity of java byte array.
            for(int iter = 0; iter < scalars.size(); iter+= G1_iteration_batch_size){
                long start = System.currentTimeMillis();
                ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );
                //we only send a single base value to CUDA for computing the windowTable per executor,
                //this is acutually faster than spark broadcast the windowTable.
                ArrayList<BigInteger> three_values = baseElement.BN254G1ToBigInteger();
                outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
                byte[] baseByteArrayXYZ = outputStream.toByteArray();
        
                final int outerc = (scalarSize + windowSize - 1) / windowSize;
                
                ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        
                for (int i = iter; i < Integer.min(iter + G1_iteration_batch_size, scalars.size()); i++) {
                    bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i).toBigInteger()));
                }
                byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
        
                long finish = System.currentTimeMillis();
                long timeElapsed = finish - start;
                //System.out.println("iteration=" + iter + " data transfer preparation time elapsed: " + timeElapsed + " ms");
                byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, Integer.min(G1_iteration_batch_size, scalars.size() - iter), scalarSize, baseByteArrayXYZ, bigScalarByteArray, 1, 0);
                
                start = System.currentTimeMillis();
                int size_of_bigint_cpp_side = 64;
                
                
                for(int i = 0; i < Integer.min(G1_iteration_batch_size, scalars.size() - iter); i++){
                    byte[] converted_back_X = Arrays.copyOfRange(resultByteArray, 3*i*size_of_bigint_cpp_side, (3*i+1)*size_of_bigint_cpp_side);
                    byte[] converted_back_Y = Arrays.copyOfRange(resultByteArray, (3*i +1)*size_of_bigint_cpp_side, (3*i+2)*size_of_bigint_cpp_side);
                    byte[] converted_back_Z = Arrays.copyOfRange(resultByteArray, (3*i +2)*size_of_bigint_cpp_side, (3*i+3)*size_of_bigint_cpp_side);
        
        
                    BigInteger bi_X = new BigInteger(converted_back_X);
                    BigInteger bi_Y = new BigInteger(converted_back_Y);
                    BigInteger bi_Z = new BigInteger(converted_back_Z);
        
        
                    T temp = baseElement.zero();
                    temp.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
                    jni_res.add(temp);
                }
        
                finish = System.currentTimeMillis();
                timeElapsed = finish - start;
                System.out.println("iteration=" + iter + " data receive transformation time elapsed: " + timeElapsed + " ms");
            }

            return jni_res;
        }else{
            //BN254 G2
            final List<T> jni_res = new ArrayList<>(scalars.size());
            int G2_iteration_batch_size = (1 << 22);

            //because when size is too large, we will exceed the capacity of java byte array.
            for(int iter = 0; iter < scalars.size(); iter += G2_iteration_batch_size){
                ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );
                ArrayList<BigInteger> six_values = baseElement.BN254G2ToBigInteger();
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
        
        
                byte[] baseByteArrayXYZ = outputStream.toByteArray();
        
                final int outerc = (scalarSize + windowSize - 1) / windowSize;
                
                ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        
                for (int i = iter ; i < Integer.min(iter + G2_iteration_batch_size, scalars.size()); i++) {
                    bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i).toBigInteger()));
                }
                byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
        
                byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, Integer.min(G2_iteration_batch_size, scalars.size() - iter), scalarSize, baseByteArrayXYZ, bigScalarByteArray, 2, 0);
                
                int size_of_bigint_cpp_side = 64;            
        
                for(int i = 0; i <  Integer.min(G2_iteration_batch_size, scalars.size() - iter); i++){
                    byte[] converted_back_Xa = Arrays.copyOfRange(resultByteArray, 6*i*size_of_bigint_cpp_side, (6*i+1)*size_of_bigint_cpp_side);
                    byte[] converted_back_Xb = Arrays.copyOfRange(resultByteArray, (6*i +1)*size_of_bigint_cpp_side, (6*i+2)*size_of_bigint_cpp_side);
                    byte[] converted_back_Ya = Arrays.copyOfRange(resultByteArray, (6*i +2)*size_of_bigint_cpp_side, (6*i+3)*size_of_bigint_cpp_side);
                    byte[] converted_back_Yb = Arrays.copyOfRange(resultByteArray, (6*i + 3)*size_of_bigint_cpp_side, (6*i+4)*size_of_bigint_cpp_side);
                    byte[] converted_back_Za = Arrays.copyOfRange(resultByteArray, (6*i +4)*size_of_bigint_cpp_side, (6*i+5)*size_of_bigint_cpp_side);
                    byte[] converted_back_Zb = Arrays.copyOfRange(resultByteArray, (6*i +5)*size_of_bigint_cpp_side, (6*i+6)*size_of_bigint_cpp_side);
        
        
                    BigInteger bi_Xa = new BigInteger(converted_back_Xa);
                    BigInteger bi_Ya = new BigInteger(converted_back_Ya);
                    BigInteger bi_Za = new BigInteger(converted_back_Za);
                    BigInteger bi_Xb = new BigInteger(converted_back_Xb);
                    BigInteger bi_Yb = new BigInteger(converted_back_Yb);
                    BigInteger bi_Zb = new BigInteger(converted_back_Zb);
        
                    T temp = baseElement.zero();
                    temp.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);
                    jni_res.add(temp);
                }
            }

            

            return jni_res;
        }
        
    

    }



    //this batchMSMPartition is for distributed zkSNARK, each executor could work on multiple partition and each partition map this method.
    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<Tuple2<Long, T>> batchMSMPartition(
            final int scalarSize,
            final int windowSize,
            final int out_size, final int in_size,
            final T baseElement,
            final List<Tuple2< Long, FieldT>> scalars, 
            final int taskID) throws Exception {
        

            if(baseElement.getClass().getName() == "algebra.curves.barreto_naehrig.bn254a.BN254aG1"){
                long start = System.currentTimeMillis();

                int G1_iteration_batch_size = (1 << 23);
                final List<Tuple2<Long, T>> jni_res = new ArrayList<>(scalars.size());
                for(int iter = 0; iter < scalars.size(); iter+= G1_iteration_batch_size){
                
                    ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );
                    //we only send a single base value to CUDA for computing the windowTable per executor,
                    //this is acutually faster than spark broadcast the windowTable.
                    ArrayList<BigInteger> three_values = baseElement.BN254G1ToBigInteger();
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
        
                    byte[] baseByteArrayXYZ = outputStream.toByteArray();
        
                    final int outerc = (scalarSize + windowSize - 1) / windowSize;
                    
                    ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
        
                    for (int i = iter ; i < Integer.min(iter + G1_iteration_batch_size, scalars.size()); i++) {
                        bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i)._2.toBigInteger()));
                    }
                    byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
        
                    long finish = System.currentTimeMillis();
                    long timeElapsed = finish - start;
                    //System.out.println("iteration=" + iter + "  data transfer preparation time elapsed: " + timeElapsed + " ms");
                    byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, Integer.min( G1_iteration_batch_size, scalars.size() - iter) , scalarSize, baseByteArrayXYZ, bigScalarByteArray, 1, (int)taskID);
                    
                    start = System.currentTimeMillis();
                    int size_of_bigint_cpp_side = 64;
                    
        
                    for(int i = 0; i < Integer.min( G1_iteration_batch_size, scalars.size() - iter); i++){
                        byte[] converted_back_X = Arrays.copyOfRange(resultByteArray, 3*i*size_of_bigint_cpp_side, (3*i+1)*size_of_bigint_cpp_side);
                        byte[] converted_back_Y = Arrays.copyOfRange(resultByteArray, (3*i +1)*size_of_bigint_cpp_side, (3*i+2)*size_of_bigint_cpp_side);
                        byte[] converted_back_Z = Arrays.copyOfRange(resultByteArray, (3*i +2)*size_of_bigint_cpp_side, (3*i+3)*size_of_bigint_cpp_side);
            
        
                        BigInteger bi_X = new BigInteger(converted_back_X);
                        BigInteger bi_Y = new BigInteger(converted_back_Y);
                        BigInteger bi_Z = new BigInteger(converted_back_Z);
        
                        T temp = baseElement.zero();
                        temp.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
                        jni_res.add(new Tuple2<>(scalars.get(i + iter)._1, temp));
                        //System.out.println("CUDA FixedBaseMSM output=" +temp.toString());
                    }
        
                    finish = System.currentTimeMillis();
                    timeElapsed = finish - start;
                    //System.out.println("iteration=" + iter + "  data receive transformation time elapsed: " + timeElapsed + " ms");
                }
                
                return jni_res;
            }else{
                ///BN254 G2 curve 
                final List<Tuple2<Long, T>> jni_res = new ArrayList<>(scalars.size());
                int G2_iteration_batch_size = (1 << 22);
                //because when size is too large, we will exceed the capacity of java byte array.
                for(int iter = 0; iter < scalars.size(); iter += G2_iteration_batch_size){
                    ByteArrayOutputStream outputStream = new ByteArrayOutputStream( );
                    ArrayList<BigInteger> six_values = baseElement.BN254G2ToBigInteger();
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                    outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
            
            
                    byte[] baseByteArrayXYZ = outputStream.toByteArray();
            
                    final int outerc = (scalarSize + windowSize - 1) / windowSize;
                    
                    ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );
                    for (int i = iter ; i < Integer.min(iter + G2_iteration_batch_size, scalars.size()); i++) {
                        bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i)._2.toBigInteger()));
                    }

                    byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
            
                    byte[] resultByteArray = batchMSMNativeHelper(outerc, windowSize, out_size, in_size, Integer.min(G2_iteration_batch_size, scalars.size() - iter), scalarSize, baseByteArrayXYZ, bigScalarByteArray, 2, (int)taskID);
                    
                    int size_of_bigint_cpp_side = 64;            
            
                    for(int i = 0; i < Integer.min(G2_iteration_batch_size, scalars.size() - iter); i++){
                        byte[] converted_back_Xa = Arrays.copyOfRange(resultByteArray, 6*i*size_of_bigint_cpp_side, (6*i+1)*size_of_bigint_cpp_side);
                        byte[] converted_back_Xb = Arrays.copyOfRange(resultByteArray, (6*i +1)*size_of_bigint_cpp_side, (6*i+2)*size_of_bigint_cpp_side);
                        byte[] converted_back_Ya = Arrays.copyOfRange(resultByteArray, (6*i +2)*size_of_bigint_cpp_side, (6*i+3)*size_of_bigint_cpp_side);
                        byte[] converted_back_Yb = Arrays.copyOfRange(resultByteArray, (6*i + 3)*size_of_bigint_cpp_side, (6*i+4)*size_of_bigint_cpp_side);
                        byte[] converted_back_Za = Arrays.copyOfRange(resultByteArray, (6*i +4)*size_of_bigint_cpp_side, (6*i+5)*size_of_bigint_cpp_side);
                        byte[] converted_back_Zb = Arrays.copyOfRange(resultByteArray, (6*i +5)*size_of_bigint_cpp_side, (6*i+6)*size_of_bigint_cpp_side);
            
            
                        BigInteger bi_Xa = new BigInteger(converted_back_Xa);
                        BigInteger bi_Ya = new BigInteger(converted_back_Ya);
                        BigInteger bi_Za = new BigInteger(converted_back_Za);
                        BigInteger bi_Xb = new BigInteger(converted_back_Xb);
                        BigInteger bi_Yb = new BigInteger(converted_back_Yb);
                        BigInteger bi_Zb = new BigInteger(converted_back_Zb);
            
                        T temp = baseElement.zero();
                        temp.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);
                        jni_res.add(new Tuple2<>(scalars.get(i + iter)._1, temp));
                    }
                }
                
                return jni_res;
        }
        
    }


    public static <GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, GroupT>
    distributedBatchMSM(
            final int scalarSize,
            final int windowSize,
            final int out_size, final int in_size,
            GroupT baseElement,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {
        final Broadcast<GroupT> baseBroadcast = sc.broadcast(baseElement);
        JavaPairRDD<Long, GroupT> true_result = scalars.mapPartitionsToPair(
            partition ->  {
                TaskContext tc = TaskContext.get();
                long taskID = tc.taskAttemptId();
                List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                long start2 = System.currentTimeMillis();
                List<Tuple2<Long, GroupT>> res =  batchMSMPartition(scalarSize, windowSize, out_size, in_size, baseBroadcast.getValue(), scalar_partition, (int)taskID);
                long finish2 = System.currentTimeMillis();
                long timeElapsed2 = finish2 - start2;
                //System.out.println("partition processing time elapsed: " + timeElapsed2 + " ms");
                return res.iterator();
            }
        );

        return true_result;

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
            final byte[] bigScalarsArrays,
            final int taskID);

    public static <G1T extends AbstractGroup<G1T>,
            G2T extends AbstractGroup<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> List<Tuple2<G1T, G2T>> doubleBatchMSM(
            final int out_size1, final int in_size1,
            final int out_size2, final int in_size2,
            final int scalarSize1,
            final int windowSize1,
            final G1T baseG1,
            final int scalarSize2,
            final int windowSize2,
            final G2T baseG2,
            final List<FieldT> scalars) throws Exception {


        // System.out.println("doubleBatchMSM info:");
        // System.out.println("batchMSM size :" + scalars.size());
        // System.out.println("scalarSize len :" + scalarSize1 + " " + scalarSize2 );
        // System.out.println("windowSize len :" + windowSize1 + " " + windowSize2);
        // System.out.println("scalars type : " + scalars.get(0).getClass().getName());

        long start = System.currentTimeMillis();

        final List<Tuple2<G1T, G2T>> jni_res = new ArrayList<>(scalars.size());

        int double_batch_iteration_batch_size = (1 << 21);



        for(int iter = 0; iter < scalars.size(); iter += double_batch_iteration_batch_size){
            ByteArrayOutputStream G1outputStream = new ByteArrayOutputStream( );

            ByteArrayOutputStream G2outputStream = new ByteArrayOutputStream( );
            ArrayList<BigInteger> three_values = baseG1.BN254G1ToBigInteger();
            G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
            G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
            G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));

            byte[] baseByteArrayXYZ = G1outputStream.toByteArray();

            ArrayList<BigInteger> six_values = baseG2.BN254G2ToBigInteger();
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
            G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));

            byte[] baseByteArrayXYZABC = G2outputStream.toByteArray();

            final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
            final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;
            ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );

            for (int i = iter ; i < Integer.min(iter + double_batch_iteration_batch_size, scalars.size()); i++) {
                bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i).toBigInteger()));
            }

            long finish = System.currentTimeMillis();
            byte[] bigScalarByteArray =  bigScalarStream.toByteArray();

            long timeElapsed = finish - start;
            System.out.println("JAVA side prepare data time elapsed: " + timeElapsed + " ms");

            byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2,
                                            out_size1, in_size1, out_size2, in_size2, Integer.min(double_batch_iteration_batch_size, scalars.size() - iter),
                                            baseByteArrayXYZ, baseByteArrayXYZABC, bigScalarByteArray, 0);



            start = System.currentTimeMillis();
            int size_of_bigint_cpp_side = 64;

            // because each G1 value takes up 3 BigIntegers, and each G2 takes up 6 BigIntegers.
            for(int i = 0; i < Integer.min(double_batch_iteration_batch_size, scalars.size() - iter); i++){

                byte[] converted_back_X = Arrays.copyOfRange(resultByteArray, 9*i*size_of_bigint_cpp_side, (9*i+1)*size_of_bigint_cpp_side);
                byte[] converted_back_Y = Arrays.copyOfRange(resultByteArray, (9*i +1)*size_of_bigint_cpp_side, (9*i+2)*size_of_bigint_cpp_side);
                byte[] converted_back_Z = Arrays.copyOfRange(resultByteArray, (9*i +2)*size_of_bigint_cpp_side, (9*i+3)*size_of_bigint_cpp_side);

                BigInteger bi_X = new BigInteger(converted_back_X);
                BigInteger bi_Y = new BigInteger(converted_back_Y);
                BigInteger bi_Z = new BigInteger(converted_back_Z);
                //System.out.println("G1 X,Y,Z=" + bi_X + " " + bi_Y + " " + bi_Z);
                G1T temp1 = baseG1.zero();
                temp1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);

                byte[] converted_back_Xa = Arrays.copyOfRange(resultByteArray, (9*i +3)*size_of_bigint_cpp_side, (9*i+4)*size_of_bigint_cpp_side);
                
                byte[] converted_back_Xb = Arrays.copyOfRange(resultByteArray, (9*i +4)*size_of_bigint_cpp_side, (9*i+5)*size_of_bigint_cpp_side);
                byte[] converted_back_Ya = Arrays.copyOfRange(resultByteArray, (9*i +5)*size_of_bigint_cpp_side, (9*i+6)*size_of_bigint_cpp_side);
                byte[] converted_back_Yb = Arrays.copyOfRange(resultByteArray, (9*i +6)*size_of_bigint_cpp_side, (9*i+7)*size_of_bigint_cpp_side);
                byte[] converted_back_Za = Arrays.copyOfRange(resultByteArray, (9*i +7)*size_of_bigint_cpp_side, (9*i+8)*size_of_bigint_cpp_side);
                byte[] converted_back_Zb = Arrays.copyOfRange(resultByteArray, (9*i +8)*size_of_bigint_cpp_side, (9*i+9)*size_of_bigint_cpp_side);

                BigInteger bi_Xa = new BigInteger(converted_back_Xa);
                BigInteger bi_Ya = new BigInteger(converted_back_Ya);
                BigInteger bi_Za = new BigInteger(converted_back_Za);
                BigInteger bi_Xb = new BigInteger(converted_back_Xb);
                BigInteger bi_Yb = new BigInteger(converted_back_Yb);
                BigInteger bi_Zb = new BigInteger(converted_back_Zb);

                G2T temp2 = baseG2.zero();
                temp2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);

                jni_res.add(new Tuple2<>(temp1, temp2));
            }

            finish = System.currentTimeMillis();
            timeElapsed = finish - start;
            System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
        }

        

        return jni_res;
    }

    public static <G1T extends AbstractGroup<G1T>,
        G2T extends AbstractGroup<G2T>,
        FieldT extends AbstractFieldElementExpanded<FieldT>> List<Tuple2< Long, Tuple2<G1T, G2T>>> doubleBatchMSMPartition(
        final int out_size1, final int in_size1,
        final int out_size2, final int in_size2,
        final int scalarSize1,
        final int windowSize1,
        final G1T baseG1,
        final int scalarSize2,
        final int windowSize2,
        final G2T baseG2,
        final List<Tuple2<Long, FieldT>> scalars, 
        final int taskID) throws Exception{

            int double_batch_iteration_batch_size = (1 << 21);
            final List<Tuple2< Long, Tuple2<G1T, G2T>>> jni_res = new ArrayList<>(scalars.size());

            for(int iter = 0; iter < scalars.size(); iter += double_batch_iteration_batch_size){
                
                long start = System.currentTimeMillis();
                ByteArrayOutputStream G1outputStream = new ByteArrayOutputStream( );
        
                ByteArrayOutputStream G2outputStream = new ByteArrayOutputStream( );
        
                ArrayList<BigInteger> three_values = baseG1.BN254G1ToBigInteger();
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(0)));
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(1)));
                G1outputStream.write(bigIntegerToByteArrayHelperCGBN(three_values.get(2)));
        
                byte[] baseByteArrayXYZ = G1outputStream.toByteArray();
        
                ArrayList<BigInteger> six_values = baseG2.BN254G2ToBigInteger();
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(0)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(1)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(2)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(3)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(4)));
                G2outputStream.write(bigIntegerToByteArrayHelperCGBN(six_values.get(5)));
        
                byte[] baseByteArrayXYZABC = G2outputStream.toByteArray();
        
                final int outerc1 = (scalarSize1 + windowSize1 - 1) / windowSize1;
                final int outerc2 = (scalarSize2 + windowSize2 - 1) / windowSize2;
                ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );

                for (int i = iter ; i < Integer.min(iter + double_batch_iteration_batch_size, scalars.size()); i++) {
                    bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalars.get(i)._2.toBigInteger()));
                }
                long finish = System.currentTimeMillis();
                byte[] bigScalarByteArray =  bigScalarStream.toByteArray();
        
                long timeElapsed = finish - start;
                //System.out.println("JAVA side prepare data time elapsed: " + timeElapsed + " ms");
        
                byte[] resultByteArray = doubleBatchMSMNativeHelper(outerc1, windowSize1, outerc2, windowSize2,
                                                out_size1, in_size1, out_size2, in_size2, Integer.min(double_batch_iteration_batch_size, scalars.size() - iter),
                                                baseByteArrayXYZ, baseByteArrayXYZABC, bigScalarByteArray, taskID);
        
        
        
                start = System.currentTimeMillis();
                int size_of_bigint_cpp_side = 64;
        
                // because each G1 value takes up 3 BigIntegers, and each G2 takes up 6 BigIntegers.
                for(int i = 0; i < Integer.min(double_batch_iteration_batch_size, scalars.size() - iter); i++){
                    byte[] converted_back_X = Arrays.copyOfRange(resultByteArray, 9*i*size_of_bigint_cpp_side, (9*i+1)*size_of_bigint_cpp_side);
                    byte[] converted_back_Y = Arrays.copyOfRange(resultByteArray, (9*i +1)*size_of_bigint_cpp_side, (9*i+2)*size_of_bigint_cpp_side);
                    byte[] converted_back_Z = Arrays.copyOfRange(resultByteArray, (9*i +2)*size_of_bigint_cpp_side, (9*i+3)*size_of_bigint_cpp_side);
        
                    BigInteger bi_X = new BigInteger(converted_back_X);
                    BigInteger bi_Y = new BigInteger(converted_back_Y);
                    BigInteger bi_Z = new BigInteger(converted_back_Z);
                    //System.out.println("G1 X,Y,Z=" + bi_X + " " + bi_Y + " " + bi_Z);
                    G1T temp1 = baseG1.zero();
                    temp1.setBigIntegerBN254G1(bi_X, bi_Y, bi_Z);
        
                    byte[] converted_back_Xa = Arrays.copyOfRange(resultByteArray, (9*i +3)*size_of_bigint_cpp_side, (9*i+4)*size_of_bigint_cpp_side);
                    
                    byte[] converted_back_Xb = Arrays.copyOfRange(resultByteArray, (9*i +4)*size_of_bigint_cpp_side, (9*i+5)*size_of_bigint_cpp_side);
                    byte[] converted_back_Ya = Arrays.copyOfRange(resultByteArray, (9*i +5)*size_of_bigint_cpp_side, (9*i+6)*size_of_bigint_cpp_side);
                    byte[] converted_back_Yb = Arrays.copyOfRange(resultByteArray, (9*i +6)*size_of_bigint_cpp_side, (9*i+7)*size_of_bigint_cpp_side);
                    byte[] converted_back_Za = Arrays.copyOfRange(resultByteArray, (9*i +7)*size_of_bigint_cpp_side, (9*i+8)*size_of_bigint_cpp_side);
                    byte[] converted_back_Zb = Arrays.copyOfRange(resultByteArray, (9*i +8)*size_of_bigint_cpp_side, (9*i+9)*size_of_bigint_cpp_side);
        
                    BigInteger bi_Xa = new BigInteger(converted_back_Xa);
                    BigInteger bi_Ya = new BigInteger(converted_back_Ya);
                    BigInteger bi_Za = new BigInteger(converted_back_Za);
                    BigInteger bi_Xb = new BigInteger(converted_back_Xb);
                    BigInteger bi_Yb = new BigInteger(converted_back_Yb);
                    BigInteger bi_Zb = new BigInteger(converted_back_Zb);
        
                    G2T temp2 = baseG2.zero();
                    temp2.setBigIntegerBN254G2(bi_Xa, bi_Xb, bi_Ya, bi_Yb, bi_Za, bi_Zb);            
                    jni_res.add(new Tuple2<>(scalars.get(i + iter)._1, new Tuple2<>(temp1, temp2)));
                }
        
                finish = System.currentTimeMillis();
                timeElapsed = finish - start;
                //System.out.println("data receive transformation time elapsed: " + timeElapsed + " ms");
            }
            
    
    
    
            return jni_res;

    }

    public static <G1T extends AbstractG1<G1T>,
            G2T extends AbstractG2<G2T>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, Tuple2<G1T, G2T>>
    distributedDoubleBatchMSM(
            final int out_size1, final int in_size1,
            final int out_size2, final int in_size2,
            final int scalarSize1,
            final int windowSize1,
            final G1T baseG1,
            final int scalarSize2,
            final int windowSize2,
            final G2T baseG2,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {

        final Broadcast<G1T> baseBroadcast1 = sc.broadcast(baseG1);
        final Broadcast<G2T> baseBroadcast2 = sc.broadcast(baseG2);
        return scalars.mapPartitionsToPair(
            partition ->  {
                    TaskContext tc = TaskContext.get();
                    long taskID = tc.taskAttemptId();
                    List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                    List<Tuple2<Long, Tuple2<G1T, G2T>>> res =  doubleBatchMSMPartition(out_size1, in_size1, out_size2, in_size2,
                                                                scalarSize1, windowSize1, baseBroadcast1.getValue(), 
                                                                scalarSize2, windowSize2, baseBroadcast2.getValue(), scalar_partition, (int)taskID );
                    return res.iterator();
            }
        );

    }




    
    public static native <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    byte[] fieldBatchMSMNativeHelper(
            final byte[] bigScalarsArrays, int batch_size, int taskID);

    //this batchMSMPartition is for distributed zkSNARK, each executor could work on multiple partition and each partition map this method.
    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<Tuple2<Long, FieldT>> batchFieldMSMPartition(
            final FieldT baseElement,
            final List<Tuple2< Long, FieldT>> scalars, 
            final int taskID) throws Exception {
        
        final List<Tuple2<Long, FieldT>> jni_res = new ArrayList<>(scalars.size());
        ByteArrayOutputStream bigScalarStream = new ByteArrayOutputStream( );

        for (Tuple2<Long, FieldT> scalar : scalars) {
            bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(scalar._2.toBigInteger()));
        }

        bigScalarStream.write(bigIntegerToByteArrayHelperCGBN(baseElement.toBigInteger()));

        byte[] bigScalarByteArray =  bigScalarStream.toByteArray();

        byte[] resultByteArray = fieldBatchMSMNativeHelper( bigScalarByteArray, scalars.size(), taskID);

        int size_of_bigint_cpp_side = 64;
        for(int i = 0; i < scalars.size(); i++){
            byte[] converted_back_X = Arrays.copyOfRange(resultByteArray, i*size_of_bigint_cpp_side, (i+1)*size_of_bigint_cpp_side);
            BigInteger bi_X = new BigInteger(converted_back_X);
            FieldT tmp = baseElement.zero();
            tmp.setBigInteger(bi_X);
            jni_res.add(new Tuple2<>(scalars.get(i)._1, tmp));
        }
        return jni_res;

    }


    public static <T extends AbstractGroup<T>, FieldT extends AbstractFieldElementExpanded<FieldT>>
    List<Tuple2<Long, FieldT>> batchFilterFieldMSMPartition(
            FieldT inverseDelta, 
            FieldT inverseGamma,
            final int numInputs,
            final List<Tuple2< Long, FieldT>> scalars, 
            int type, //0 stands for gamma, 1 stands for delta.
            final int taskID) throws Exception {
        
        final List<Tuple2<Long, FieldT>> jni_res = new ArrayList<>(scalars.size());
        ByteArrayOutputStream gammaABCStream = new ByteArrayOutputStream( );
        ByteArrayOutputStream deltaABCStream = new ByteArrayOutputStream( );
        int counter = 0;
        for (Tuple2<Long, FieldT> scalar : scalars) {
            if(scalar._1 < numInputs){
                if(type == 0){
                    gammaABCStream.write(bigIntegerToByteArrayHelperCGBN(scalar._2.toBigInteger()));
                    counter += 1;
                }
            }else{
                if(type == 1){
                    deltaABCStream.write(bigIntegerToByteArrayHelperCGBN(scalar._2.toBigInteger()));
                    counter += 1;
                }
            }
        }


        if(type == 0){
            gammaABCStream.write(bigIntegerToByteArrayHelperCGBN(inverseGamma.toBigInteger()));
            byte[] gammaABCByteArray =  gammaABCStream.toByteArray();
            byte[] gammaResultByteArray = fieldBatchMSMNativeHelper( gammaABCByteArray, counter, taskID);
            int size_of_bigint_cpp_side = 64;
            int counterGamma = 0;
            for(int i = 0; i < scalars.size(); i++){
                if(scalars.get(i)._1 >= numInputs){
                    continue;
                }
                byte[] converted_back_X = Arrays.copyOfRange(gammaResultByteArray, counterGamma*size_of_bigint_cpp_side, (counterGamma+1)*size_of_bigint_cpp_side);
                BigInteger bi_X = new BigInteger(converted_back_X);
                FieldT tmp = inverseGamma.zero();
                tmp.setBigInteger(bi_X);
                jni_res.add(new Tuple2<>(scalars.get(i)._1, tmp));
                counterGamma++;
            }
        }else{
            deltaABCStream.write(bigIntegerToByteArrayHelperCGBN(inverseDelta.toBigInteger()));
            byte[] deltaABCByteArray =  deltaABCStream.toByteArray();
            byte[] deltaResultByteArray = fieldBatchMSMNativeHelper( deltaABCByteArray, counter, taskID);

            int size_of_bigint_cpp_side = 64;
            int counterDelta = 0;
            for(int i = 0; i < scalars.size(); i++){
                if(scalars.get(i)._1 < numInputs){
                    continue;
                }
                byte[] converted_back_X = Arrays.copyOfRange(deltaResultByteArray, counterDelta*size_of_bigint_cpp_side, (counterDelta+1)*size_of_bigint_cpp_side);
                BigInteger bi_X = new BigInteger(converted_back_X);
                FieldT tmp = inverseGamma.zero();
                tmp.setBigInteger(bi_X);
                jni_res.add(new Tuple2<>(scalars.get(i)._1, tmp));
                counterDelta++;
            }
        }


        return jni_res;

    }

    public static <GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, FieldT>
    distributedFilterFieldBatchMSM(
            FieldT inverseDelta, 
            FieldT inverseGamma,
            final int numInputs,
            final int type, // 0 is gamma, 1 is delta.
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {
        final Broadcast<FieldT> inverseDeltaBroadcast = sc.broadcast(inverseDelta);
        final Broadcast<FieldT> inverseGammaBroadcast = sc.broadcast(inverseGamma);


        JavaPairRDD<Long, FieldT> true_result = scalars.mapPartitionsToPair(
            partition ->  {
                TaskContext tc = TaskContext.get();
                long taskID = tc.taskAttemptId();
                List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                List<Tuple2<Long, FieldT>> res = batchFilterFieldMSMPartition( inverseDeltaBroadcast.getValue(), inverseGammaBroadcast.getValue(),  numInputs,  scalar_partition, type,(int)taskID);
                return res.iterator();
            }
        );
        return true_result;

    }

    public static <GroupT extends AbstractGroup<GroupT>,
            FieldT extends AbstractFieldElementExpanded<FieldT>> JavaPairRDD<Long, FieldT>
    distributedFieldBatchMSM(
            FieldT baseElement,
            final JavaPairRDD<Long, FieldT> scalars,
            final JavaSparkContext sc) {
        final Broadcast<FieldT> baseBroadcast = sc.broadcast(baseElement);


        JavaPairRDD<Long, FieldT> true_result = scalars.mapPartitionsToPair(
            partition ->  {
                TaskContext tc = TaskContext.get();
                long taskID = tc.taskAttemptId();
                List<Tuple2<Long, FieldT>> scalar_partition = IteratorUtils.toList(partition);
                List<Tuple2<Long, FieldT>> res = batchFieldMSMPartition( baseBroadcast.getValue(), scalar_partition, (int)taskID);
                return res.iterator();
            }
        );
        return true_result;

    }
}


