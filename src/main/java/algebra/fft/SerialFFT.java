package algebra.fft;

import algebra.fields.AbstractFieldElementExpanded;
import common.MathUtils;

import java.io.Serializable;
import java.util.List;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Arrays;

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
    /**
     * Compute the FFT, over the domain S, of the vector input, and stores the result in input.
     */
    public void radix2FFT(final List<FieldT> input) {
        assert (input.size() == domainSize);
        ArrayList<byte[]> inputByteArray = new ArrayList<byte[]>();
        for(FieldT f : input){
            inputByteArray.add(bigIntegerToByteArrayHelper(f.toBigInteger()));
        }
        byte[] omegaBytes = bigIntegerToByteArrayHelper(omega.toBigInteger());
        System.out.println("FFT input type : " + input.get(0).getClass().getName());

        FFTAuxiliary.serialRadix2FFT(input, omega);



        // System.out.println("radix2FFT java side is about to launch");
        // byte[] resultByteArray = FFTAuxiliary.serialRadix2FFTNativeHelper(inputByteArray, omegaBytes);
        // int size_of_bigint_cpp_side = 64;
        // BigInteger modulus = new BigInteger("1532495540865888858358347027150309183618765510462668801");
        // for(int i = 0; i < input.size(); i++){
        //     byte[] slice = Arrays.copyOfRange(resultByteArray, i*size_of_bigint_cpp_side, (i+1)*size_of_bigint_cpp_side);//in cpp side, BigInt is 32 bytes.

        //     byte[] converted_back = new byte[64];
        //     for(int j = 63; j >= 3; j-=4){
        //         converted_back[j] = slice[j - 3];
        //         converted_back[j-1] = slice[j - 2];
        //         converted_back[j-2] = slice[j - 1];
        //         converted_back[j-3] = slice[j ];
        //     }

        //     BigInteger bi = new BigInteger(converted_back);
        //     BigInteger output = bi.mod(modulus);
        //     FieldT temp = input.get(0).zero();
        //     temp.setBigInteger(output);
        //     input.set(i, temp);

        // }
    }

    /**
     * Compute the inverse FFT, over the domain S, of the vector input, and stores the result in
     * input.
     */
    public  void radix2InverseFFT(final List<FieldT> input) {
        //TODO lianke implement JNI code for it
        assert (input.size() == domainSize);
        System.out.println("radix2InverseFFT input side " + input.size());
        FFTAuxiliary.serialRadix2FFT(input, omega.inverse());

        final FieldT constant = input.get(0).construct(domainSize).inverse();
        for (int i = 0; i < domainSize; ++i) {
            input.set(i, input.get(i).mul(constant));
        }
    }

    /**
     * Compute the FFT, over the domain g*S, of the vector input, and stores the result in input.
     */
    public void radix2CosetFFT(final List<FieldT> input, final FieldT g) {
        System.out.println("radix2CosetFFT input side " + input.size());

        //TODO lianke implement JNI code for it
        FFTAuxiliary.multiplyByCoset(input, g);
        this.radix2FFT(input);
    }

    /**
     * Compute the inverse FFT, over the domain g*S, of the vector input, and stores the result in
     * input.
     */
    public void radix2CosetInverseFFT(final List<FieldT> input, final FieldT g) {
        //TODO lianke implement JNI code for it
        System.out.println("radix2CosetInverseFFT input side " + input.size());
        this.radix2InverseFFT(input);
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

