package scs;


public class TestProblem {
    public static void main(String [] args) {
        int m = 50; // rows
        int n = 30; // cols

        AMatrix A = AMatrix.generateRandomMatrix(m, n);
        double[] b = Utils.generateRandomDoubleArray(m);
        double[] c = Utils.generateRandomDoubleArray(n);
        Cone k = new Cone();
        k.setL(m); // random LP

        Data d = new Data(A, b, c);
        Params p = new Params();
        ConeProgram cp;
        Solution sol;

        cp = new ConeProgram(d, k, p, new DirectSolver());
        sol = cp.solve();
        System.out.println("c'x = " + Utils.ip(sol.getX(), c));
        System.out.println("b'y = " + Utils.ip(sol.getY(), b));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(A, b, sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(A, c, sol));

        cp = new ConeProgram(d, k, p, new IndirectSolver());
        sol = cp.solve();
        System.out.println("c'x = " + Utils.ip(sol.getX(), c));
        System.out.println("b'y = " + Utils.ip(sol.getY(), b));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(A, b, sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(A, c, sol));
    }
}
