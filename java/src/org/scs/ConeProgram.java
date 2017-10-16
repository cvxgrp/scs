package org.scs;

public class ScsConeProgram {
    protected ScsData d;
    protected ScsCone k;
    protected ScsSettings p;
    protected IScsConeSolver solver;

    ScsConeProgram() {}

    public ScsConeProgram(ScsData d, ScsCone k, ScsSettings p, IScsConeSolver solver) {
        this.d = d;
        this.k = k;
        this.p = p;
        this.solver = solver;
    }

    public ScsConeProgram(ScsMatrix A, double[] b, double[] c, ScsCone k, ScsSettings p, IScsConeSolver solver) {
        this(new ScsData(A, b, c), k, p, solver);
    }

    public ScsSolutionution solve() {
        ScsSolutionution sol = new ScsSolutionution();
        ScsInfo info = new ScsInfo();
        solver.solve(d, k, p, sol, info);
        sol.setScsInfo(info);
        return sol;
    }
}
