package scs;

public class ConeProgram {
    private Data d;
    private Cone k;
    private Params p;

    static {
        //System.loadLibrary("scsindir"); // Load native library at runtime
        System.loadLibrary("scsdir"); // Load native library at runtime
    }
    private static native Solution csolve(AMatrix A, double[] b, double[] c, Cone k, Params p, Solution s);

    public ConeProgram(Data d, Cone k, Params p) {
        this.d = d;
        this.k = k;
        this.p = p;
    }

    public ConeProgram(AMatrix A, double[] b, double[] c, Cone k, Params p) {
        this.d = new Data(A, b, c);
        this.k = k;
        this.p = p;
    }

    public Solution solve() {
        Solution sol = new Solution();
        csolve(d.getA(), d.getB(), d.getC(), k, p, sol);
        return sol;
    }
}
