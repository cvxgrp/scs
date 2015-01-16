package scs;

public class SolutionWithInfo {
    private final Solution sol;
    private final Info info;

    public SolutionWithInfo(Solution sol, Info info) {
        this.sol = sol;
        this.info = info;
    }

    public Solution getSol() {
        return sol;
    }

    public Info getInfo() {
        return info;
    }
}