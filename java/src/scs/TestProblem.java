package scs;


public class TestProblem {
    public static void main(String [] args) {
        int m;
        int n;
        if (args.length < 2) {
            m = 50; // rows
            n = 30; // cols
        } else {
            m = Integer.parseInt(args[0]);
            n = Integer.parseInt(args[1]);
        }

        Params p = new Params();
        RandomLinearProgram cp;
        Solution sol;
        Data d = new Data();

        cp = new RandomLinearProgram(m, n, d, p, new DirectSolver());
        System.out.println("true opt = " + cp.getOpt());
        sol = cp.solve();

        System.out.println("c'x = " + Utils.ip(sol.getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(sol.getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), sol));

        cp.setSolver(new IndirectSolver());
        sol = cp.solve();

        System.out.println("c'x = " + Utils.ip(sol.getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(sol.getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), sol));
    }
}
