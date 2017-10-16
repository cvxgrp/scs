package org.scs;

public class DirectSolver implements IScsConeSolver {
    static {
        System.loadLibrary("jscsdir"); // Load native library at runtime
    }

    private static native void csolve(ScsMatrix A, double[] b, double[] c, ScsCone k, ScsSettings p, ScsSolutionution s, ScsInfo info);
    private static native String cversion();

    private final static String VERSION = cversion();

    public void solve(ScsData d, ScsCone k, ScsSettings p, ScsSolutionution sol, ScsInfo info) {
        csolve(d.getA(), d.getB(), d.getC(), k, p, sol, info);
    }

    public String version() {
        return VERSION;
    }
}
