package org.scs;

public class ScsSolutionution {
    private ScsInfo info;  // information about the run
    private double[] x; // primal variable
    private double[] y; // dual variable
    private double[] s; // primal slack variable

    public ScsInfo getScsInfo() {
        return info;
    }
    public double[] getX() {
        return x;
    }
    public double[] getY() {
        return y;
    }
    public double[] getS() {
        return s;
    }

    public void setScsInfo(ScsInfo info) {
        this.info = info;
    }
    public void setX(double[] x) {
        this.x = x;
    }
    public void setY(double[] y) {
        this.y = y;
    }
    public void setS(double[] s) {
        this.s = s;
    }
}