package org.scs;

public interface IConeSolver {
    public void solve(Data d, Cone k, Settings p, Solution sol, Info info);
    public String version();
}
