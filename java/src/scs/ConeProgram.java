package scs;

public class ConeProgram {
    protected Data d;
    protected Cone k;
    protected Params p;
    protected IConeSolver solver;

    ConeProgram() {};

    public ConeProgram(Data d, Cone k, Params p, IConeSolver solver) {
        this.d = d;
        this.k = k;
        this.p = p;
        this.solver = solver;
    }

    public ConeProgram(AMatrix A, double[] b, double[] c, Cone k, Params p, IConeSolver solver) {
        this(new Data(A, b, c), k, p, solver);
    }

    public Solution solve() {
        return solver.solve(d, k, p);
    }
}
