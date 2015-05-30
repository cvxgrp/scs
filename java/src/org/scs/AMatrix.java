package org.scs;

import java.util.Random;
import java.util.TreeSet;

public class AMatrix {
    /**
     * A matrix, dense format, size m by n, column layout
     */
    private double[] v;  /* A values, size: m * n */
    private int m;
    private int n;

    public AMatrix(int m, int n, double[] v) {
        this.m = m;
        this.n = n;
        this.v = v;
    }

    public double[] getValues() {
        return v;
    }

    public int getNumRows() {
        return m;
    }

    public int getNumCols() {
        return n;
    }

    /* Utilities: */
    public void accumByA(double[] x, double[] y) {
        // y += Ax
        for (int j = 0; j < n; j++) { // cols
            for (int i = 0; i < m; i++) { // rows
                y[i] += v[j*m + i] * x[j];
            }
        }
    }

    public void accumByAtrans(double[] x, double[] y) {
        // y += A'x
        for (int j = 0; j < n; j++) { // cols
            for (int i = 0; i < m; i++) { // rows
                y[j] += v[j*m + i] * x[i];
            }
        }
    }
    
    public static AMatrix generateRandomMatrix(int m, int n) {
        AMatrix A = new AMatrix(m, n, new double[n * m]);
        Random rng = new Random();
        for (int j = 0; j < n * m; j++) {
            A.v[j] = 10 * (rng.nextDouble() - 0.5);
        }
        return A;
    }
}
