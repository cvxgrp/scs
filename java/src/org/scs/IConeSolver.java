package org.scs;

public interface IScsConeSolver {
    public void solve(ScsData d, ScsCone k, ScsSettings p, ScsSolutionution sol, ScsInfo info);
    public String version();
}
