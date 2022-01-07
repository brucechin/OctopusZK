package algebra.fft;

import algebra.fields.AbstractFieldElementExpanded;
import common.MathUtils;

import java.io.Serializable;
import java.util.List;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;
import scala.Tuple2;

import algebra.groups.AbstractGroup;

public class SerialFFT<FieldT extends AbstractFieldElementExpanded<FieldT>> implements
        Serializable {

    public final int domainSize;
    private final FieldT omega;
    static {
		System.loadLibrary("AlgebraFFTAuxiliary");
        System.out.println("AlgebraFFTAuxiliary loaded");
	}
    public SerialFFT(final long _domainSize, final FieldT fieldFactory) {
        assert (_domainSize > 1);
        domainSize = (int) MathUtils.lowestPowerOfTwo(_domainSize);
        omega = fieldFactory.rootOfUnity(domainSize);
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
     * Compute the FFT, over the domain S, of the vector input, and stores the result in input.
     */
    public void radix2FFT(final List<FieldT> input, final int taskID) throws Exception{
        assert (input.size() == domainSize);
        FFTAuxiliary.serialRadix2FFT(input, omega, taskID);
    }

    public void radix2FFTBatch(final List<Tuple2<Long, ArrayList<Tuple2<Long, FieldT>>>>  input, final int taskID) throws Exception{
        assert (input.size() == domainSize);
        FFTAuxiliary.serialRadix2FFTBatchPartition(input, omega, taskID);
    }


    /**
     * Compute the inverse FFT, over the domain S, of the vector input, and stores the result in
     * input.
     */
    public  void radix2InverseFFT(final List<FieldT> input, final int taskID) throws Exception{
        assert (input.size() == domainSize);
        //System.out.println("radix2InverseFFT input side " + input.size());
        FFTAuxiliary.serialRadix2FFT(input, omega.inverse(), taskID);

        //TODO lianke maybe move this to GPU
        final FieldT constant = input.get(0).construct(domainSize).inverse();
        for (int i = 0; i < domainSize; ++i) {
            input.set(i, input.get(i).mul(constant));
        }
    }

    public  void radix2InverseFFTBatch(final List<Tuple2<Long, ArrayList<Tuple2<Long, FieldT>>>> input, final int taskID) throws Exception{
        assert (input.get(0)._2.size() == domainSize);
        System.out.println("radix2InverseFFT input side " + input.size() + " domain size " + input.get(0)._2.size());
        FFTAuxiliary.serialRadix2FFTBatchPartition(input, omega.inverse(), taskID);

        //TODO lianke maybe move this to GPU
        for(int i = 0; i < input.size(); i++){
            final FieldT constant = input.get(i)._2.get(0)._2.construct(domainSize).inverse();
            for (int j = 0; j < domainSize; j++) {
                input.get(i)._2.set(j, new Tuple2<>(input.get(i)._2.get(j)._1, input.get(i)._2.get(j)._2.mul(constant)));
            }
        }

    }

    /**
     * Compute the FFT, over the domain g*S, of the vector input, and stores the result in input.
     */
    public void radix2CosetFFT(final List<FieldT> input, final FieldT g, final int taskID)throws Exception {
        //System.out.println("radix2CosetFFT input side " + input.size());

        FFTAuxiliary.multiplyByCoset(input, g);
        this.radix2FFT(input, taskID);
    }

    /**
     * Compute the inverse FFT, over the domain g*S, of the vector input, and stores the result in
     * input.
     */
    public void radix2CosetInverseFFT(final List<FieldT> input, final FieldT g, final int taskID)throws Exception {
        //System.out.println("radix2CosetInverseFFT input side " + input.size());
        this.radix2InverseFFT(input, taskID);
        FFTAuxiliary.multiplyByCoset(input, g.inverse());
    }

    /**
     * Evaluate all Lagrange polynomials.
     * <p>
     * The inputs are:
     * - an integer m
     * - an element t
     * The output is a vector (b_{0},...,b_{m-1})
     * where b_{i} is the evaluation of L_{i,S}(z) at z = t.
     */
    public List<FieldT> lagrangeCoefficients(final FieldT t) {
        return FFTAuxiliary.serialRadix2LagrangeCoefficients(t, domainSize);
    }

    /**
     * Returns the ith power of omega, omega^i.
     */
    public FieldT getDomainElement(final int i) {
        return omega.pow(i);
    }

    /**
     * Evaluate the vanishing polynomial of S at the field element t.
     */
    public FieldT computeZ(final FieldT t) {
        return t.pow(domainSize).sub(t.one());
    }

    /**
     * Add the coefficients of the vanishing polynomial of S to the coefficients of the polynomial H.
     */
    public void addPolyZ(final FieldT coefficient, final List<FieldT> H) {
        assert (H.size() == domainSize + 1);
        H.set(domainSize, H.get(domainSize).add(coefficient));
        H.set(0, H.get(0).sub(coefficient));
    }

    /**
     * Multiply by the evaluation, on a coset of S, of the inverse of the vanishing
     * polynomial of S, and stores the result in input.
     */
    public void divideByZOnCoset(final FieldT coset, final List<FieldT> input) {
        final FieldT inverseZCoset = computeZ(coset).inverse();
        for (int i = 0; i < domainSize; i++) {
            input.set(i, input.get(i).mul(inverseZCoset));
        }
    }

}

