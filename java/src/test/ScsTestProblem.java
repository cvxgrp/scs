import org.scs.*;

public class ScsTestProblem {
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
        Solution sol;
        Data d = new Data();
        IConeSolver isolver = new IndirectSolver();
        IConeSolver dsolver = new DirectSolver();
        RandomLinearProgram cp = new RandomLinearProgram(m, n, d, p, isolver);

        sol = cp.solve();
        System.out.println("solver version: " + isolver.version());
        System.out.println("true opt = " + cp.getOpt());
        System.out.println("c'x = " + Utils.ip(sol.getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(sol.getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), sol));

        cp.setSolver(dsolver); // test direct solver

        sol = cp.solve();
        System.out.println("solver version: " + dsolver.version());
        System.out.println("true opt = " + cp.getOpt());
        System.out.println("c'x = " + Utils.ip(sol.getX(), d.getC()));
        System.out.println("b'y = " + Utils.ip(sol.getY(), d.getB()));
        System.out.println("||Ax + s - b|| / (1 + ||b||)  = " + Utils.getScaledPriResidNorm(d.getA(), d.getB(), sol));
        System.out.println("||A'y + c|| / (1 + ||c||)  = " + Utils.getScaledDualResidNorm(d.getA(), d.getC(), sol));

        /* extra info */
        System.out.println("iters " + sol.getInfo().getIter());
        System.out.println("status " + sol.getInfo().getStatus());
        System.out.println("pobj " + sol.getInfo().getPobj());
        System.out.println("dobj " + sol.getInfo().getDobj());
        System.out.println("resPri " + sol.getInfo().getResPri());
        System.out.println("resDual " + sol.getInfo().getResDual());
        System.out.println("relGap " + sol.getInfo().getRelGap());
        System.out.println("setup time " + sol.getInfo().getSetupTime());
    }
}
