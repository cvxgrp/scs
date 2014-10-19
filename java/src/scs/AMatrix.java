package scs;

import java.util.Random;

public class AMatrix {
    /**
     * A matrix, column compressed sparse format, size m by n
     */
    public double[] x;  /* A values, size: NNZ A */
    public int[] i;     /* A row index, size: NNZ A */
    public int[] p;     /* A col index, size: n+1 */

    public AMatrix(double[] x, int[] i, int[] p) {
        this.x = x;
        this.i = i;
        this.p = p;
    }

    public double[] getValues() {
        return x;
    }

    public int[] getRowIdxs() {
        return i;
    }

    public int[] getColIdxs() {
        return p;
    }

    /* Utilities for random matrix generation: */
    public static AMatrix generateRandomMatrix(int m, int n) {
        return generateRandomMatrix(m, n, 0.25f);
    }

    public static AMatrix generateRandomMatrix(int m, int n, float density) {
        int nnz = (int) (n * m * density);
        int col_nnz = Math.min((int) (nnz / n), m);
        return generateRandomMatrix(m, n, nnz, col_nnz);
    }

    public static AMatrix generateRandomMatrix(int m, int n, int nnz, int col_nnz) {
        AMatrix A = new AMatrix(new double[nnz], new int[nnz], new int[n + 1]);
        Random rng = new Random();
        A.p[0] = 0;
        for (int j = 0; j < n; j++) { /* column */
            for (int r = 0; r < col_nnz; r++) { /* row index */
                int i = Math.abs(rng.nextInt()) % m; /* row */
                A.x[r + j * col_nnz] = rng.nextDouble();
                A.i[r + j * col_nnz] = i;
            }
            A.p[j + 1] = (j + 1) * col_nnz;
        }
        return A;
    }
}