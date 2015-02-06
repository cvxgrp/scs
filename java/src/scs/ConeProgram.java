package scs;

public class ConeProgram {
    protected Data d;
    protected Cone k;
    protected Settings p;
    protected IConeSolver solver;

    ConeProgram() {}

    public ConeProgram(Data d, Cone k, Settings p, IConeSolver solver) {
        this.d = d;
        this.k = k;
        this.p = p;
        this.solver = solver;
    }

    public ConeProgram(AMatrix A, double[] b, double[] c, Cone k, Settings p, IConeSolver solver) {
        this(new Data(A, b, c), k, p, solver);
    }

    public SolutionWithInfo solve() {
        Solution sol = new Solution();
        Info info = new Info();
        solver.solve(d, k, p, sol, info);
        return new SolutionWithInfo(sol, info);
    }
}
