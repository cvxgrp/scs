package scs;

public class IndirectSolver implements IConeSolver {
    static {
        System.loadLibrary("jscsindir"); // Load native library at runtime
    }
    private static native void csolve(AMatrix A, double[] b, double[] c, Cone k, Params p, Solution s);
    public Solution solve(Data d, Cone k, Params p) {
        Solution sol = new Solution();
        csolve(d.getA(), d.getB(), d.getC(), k, p, sol);
        return sol;
    }
}
