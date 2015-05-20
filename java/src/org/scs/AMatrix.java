package org.scs;

import java.util.Random;
import java.util.TreeSet;

public class AMatrix {
    /**
     * A matrix, column compressed sparse format, size m by n
     */
    private double[] v;  /* A values, size: NNZ A */
    private int[] i;     /* A row index, size: NNZ A */
    private int[] p;     /* A col ptr, size: n+1 */
    private int m;
    private int n;

    public AMatrix(int m, int n, double[] v, int[] i, int[] p) {
        this.v = v;
        this.i = i;
        this.p = p;
        this.m = m;
        this.n = n;
    }

    public double[] getValues() {
        return v;
    }

    public int[] getRowIdxs() {
        return i;
    }

    public int[] getColIdxs() {
        return p;
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
        for (int j = 0; j < n; j++) {
            for (int c = p[j]; c < p[j + 1]; c++) {
                y[i[c]] += v[c] * x[j];
            }
        }
    }

    public void accumByAtrans(double[] x, double[] y) {
        // y += A'x
        for (int j = 0; j < n; j++) {
            for (int c = p[j]; c < p[j + 1]; c++) {
                y[j] += v[c] * x[i[c]];
            }
        }
    }
    
    /* used for generating random instances */
    private static final double DEFAULT_DENSITY = 0.25;

    public static AMatrix generateRandomMatrix(int m, int n) {
        return generateRandomMatrix(m, n, DEFAULT_DENSITY);
    }

    public static AMatrix generateRandomMatrix(int m, int n, double density) {
        int col_nnz = Math.min(m, (int) (m * density));
        return generateRandomMatrix(m, n, col_nnz);
    }

    public static AMatrix generateRandomMatrix(int m, int n, int col_nnz) {
        AMatrix A = new AMatrix(m, n, new double[n * col_nnz], new int[n * col_nnz], new int[n + 1]);
        Random rng = new Random();
        A.p[0] = 0;
        for (int j = 0; j < n; j++) { /* column */
            TreeSet<Integer> rows = new TreeSet<Integer>();
            while (rows.size() < col_nnz) {
                rows.add(Math.abs(rng.nextInt()) % m); /* row */
            }
            int r = 0;
            for (Integer i : rows) {
                A.v[r + j * col_nnz] = 10 * (rng.nextDouble() - 0.5);
                A.i[r + j * col_nnz] = i;
                r++;
            }
            A.p[j + 1] = (j + 1) * col_nnz;
        }
        return A;
    }
}
