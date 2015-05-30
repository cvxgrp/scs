package org.scs;

import java.util.Random;

public class Utils {
    public static double ip(double[] a, double[] b) {
        double sum = 0;
        for (int i = 0; i < a.length; i++) {
            sum += a[i] * b[i];
        }
        return sum;
    }

    public static void addScaledArray(double[] a, double[] b, double scale) {
        for (int i = 0; i < a.length; i++) {
            b[i] += scale * a[i];
        }
    }

    public static void scaleArray(double[] a, double scale) {
        for (int i = 0; i < a.length; i++) {
            a[i] *= scale;
        }
    }

    public static double norm(double[] a) {
        double nrmSq = 0;
        for (int i = 0; i < a.length; i++) {
            nrmSq += a[i] * a[i];
        }
        return Math.sqrt(nrmSq);
    }

    public static double[] generateRandomDoubleArray(int l) {
        Random rng = new Random();
        double[] a = new double[l];
        for (int i = 0; i < l; i++) {
            a[i] = 10 * (rng.nextDouble() - 0.5);
        }
        return a;
    }

    public static double getScaledPriResidNorm(AMatrix A, double[] b, Solution sol) {
        double[] priResid = new double[A.getNumRows()];
        A.accumByA(sol.getX(), priResid);
        Utils.addScaledArray(sol.getS(), priResid, 1d);
        Utils.addScaledArray(b, priResid, -1d);
        return Utils.norm(priResid) / (1 + Utils.norm(b));
    }

    public static double getScaledDualResidNorm(AMatrix A, double[] c, Solution sol) {
        double[] dualResid = new double[A.getNumCols()];
        A.accumByAtrans(sol.getY(), dualResid);
        Utils.addScaledArray(c, dualResid, 1d);
        return Utils.norm(dualResid) / (1 + Utils.norm(c));
    }
}
