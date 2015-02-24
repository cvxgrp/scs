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

        Settings p = new Settings();
        SolutionWithInfo si;
        Data d = new Data();
        IConeSolver isolver = new IndirectSolver();
        IConeSolver dsolver = new IndirectSolver();
        RandomLinearProgram cp = new RandomLinearProgram(m, n, d, p, isolver);

        si = cp.solve();
        System.out.println("solver version: " + isolver.version());
        System.out.println("true opt = " + cp.getOpt());
        System.out.println("c'x = " + Utils.ip(si.getSol().getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(si.getSol().getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), si.getSol()));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), si.getSol()));

        cp.setSolver(dsolver); // test direct solver

        si = cp.solve();
        System.out.println("solver version: " + dsolver.version());
        System.out.println("true opt = " + cp.getOpt());
        System.out.println("c'x = " + Utils.ip(si.getSol().getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(si.getSol().getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), si.getSol()));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), si.getSol()));

        /* extra info */
        System.out.println("iters " + si.getInfo().getIter());
        System.out.println("status " + si.getInfo().getStatus());
        System.out.println("pobj " + si.getInfo().getPobj());
        System.out.println("dobj " + si.getInfo().getDobj());
        System.out.println("resPri " + si.getInfo().getResPri());
        System.out.println("resDual " + si.getInfo().getResDual());
        System.out.println("relGap " + si.getInfo().getRelGap());
        System.out.println("setup time " + si.getInfo().getSetupTime());
    }
}
