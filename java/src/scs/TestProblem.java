package scs;

import java.util.Random;

public class TestProblem {
    public static void main(String [] args) {
        int m = 50; // rows
        int n = 50; // cols

        AMatrix A = AMatrix.generateRandomMatrix(m, n, 1f);
        double[] b = generateRandomDoubleArray(m);
        double[] c = generateRandomDoubleArray(n);

        Data d = new Data(A, b, c);
        Cone k = new Cone();
        k.setF(m);
        Params p = new Params();

        ConeProgram cp = new ConeProgram(d, k, p);
        Solution sol = cp.solve();

        double sum = 0;
        for (int i=0; i < n; i++) {
            sum += sol.x[i] * c[i];
        }
        System.out.println(sum);

        sum = 0;
        for (int i=0; i < m; i++) {
            sum += sol.y[i] * b[i];
        }
        System.out.println(-sum);

    }

    private static double[] generateRandomDoubleArray(int l) {
        Random rng = new Random();
        double[] a = new double[l];
        for (int i=0; i<l; i++) {
            a[i] = rng.nextDouble();
        }
        return a;
    }
}
