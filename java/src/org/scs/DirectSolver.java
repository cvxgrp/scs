package org.scs;

public class DirectSolver implements IConeSolver {
    static {
        System.loadLibrary("jscsdir"); // Load native library at runtime
    }

    private static native void csolve(AMatrix A, double[] b, double[] c, Cone k, Settings p, Solution s, Info info);
    private static native String cversion();

    private final static String VERSION = cversion();

    public void solve(Data d, Cone k, Settings p, Solution sol, Info info) {
        csolve(d.getA(), d.getB(), d.getC(), k, p, sol, info);
    }

    public String version() {
        return VERSION;
    }
}
