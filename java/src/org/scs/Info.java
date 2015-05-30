package org.scs;

public class Info {
    private int iter; /* number of iterations taken */
    private String status; /* status string, e.g. 'Solved' */
    private int statusVal; /* status as scs_int, defined below */
    private double pobj; /* primal objective */
    private double dobj; /* dual objective */
    private double resPri; /* primal equality residual */
    private double resDual; /* dual equality residual */
    private double resInfeas; /* infeasibility cert residual */
    private double resUnbdd; /* unbounded cert residual */
    private double relGap; /* relative duality gap */
    private double setupTime; /* time taken for setup phase */
    private double solveTime; /* time taken for solve phase */

    public void setIter(int iter) {
        this.iter = iter;
    }

    public void setStatus(String status) {
        this.status = status;
    }

    public void setStatusVal(int statusVal) {
        this.statusVal = statusVal;
    }

    public void setPobj(double pobj) {
        this.pobj = pobj;
    }

    public void setDobj(double dobj) {
        this.dobj = dobj;
    }

    public void setResPri(double resPri) {
        this.resPri = resPri;
    }

    public void setResDual(double resDual) {
        this.resDual = resDual;
    }

    public void setResInfeas(double resInfeas) {
        this.resInfeas = resInfeas;
    }

    public void setResUnbdd(double resUnbdd) {
        this.resUnbdd = resUnbdd;
    }

    public void setRelGap(double relGap) {
        this.relGap = relGap;
    }

    public void setSetupTime(double setupTime) {
        this.setupTime = setupTime;
    }

    public void setSolveTime(double solveTime) {
        this.solveTime = solveTime;
    }

    public int getIter() {
        return iter;
    }

    public String getStatus() {
        return status;
    }

    public int getStatusVal() {
        return statusVal;
    }

    public double getPobj() {
        return pobj;
    }

    public double getDobj() {
        return dobj;
    }

    public double getResPri() {
        return resPri;
    }

    public double getResDual() {
        return resDual;
    }

    public double getResInfeas() {
        return resInfeas;
    }

    public double getResUnbdd() {
        return resUnbdd;
    }

    public double getRelGap() {
        return relGap;
    }

    public double getSetupTime() {
        return setupTime;
    }

    public double getSolveTime() {
        return solveTime;
    }
}
