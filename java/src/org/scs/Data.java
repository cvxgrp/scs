package org.scs;

public class ScsData {
    private ScsMatrix A;
    private double[] b;
    private double[] c;

    public ScsData() {};

    public ScsData(ScsMatrix A, double[] b, double[] c) {
        this.A = A;
        this.b = b;
        this.c = c;
    }

    // getters:
    public ScsMatrix getA() {
        return A;
    }
    public double[] getB() {
        return b;
    }
    public double[] getC() {
        return c;
    }

    // setters:
    public void setA(ScsMatrix A) {
        this.A = A;
    }
    public void setB(double[] b) {
        this.b = b;
    }
    public void setC(double[] c) {
        this.c = c;
    }
}