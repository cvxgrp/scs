package org.scs;

public class Data {
    private AMatrix A;
    private double[] b;
    private double[] c;

    public Data() {};

    public Data(AMatrix A, double[] b, double[] c) {
        this.A = A;
        this.b = b;
        this.c = c;
    }

    // getters:
    public AMatrix getA() {
        return A;
    }
    public double[] getB() {
        return b;
    }
    public double[] getC() {
        return c;
    }

    // setters:
    public void setA(AMatrix A) {
        this.A = A;
    }
    public void setB(double[] b) {
        this.b = b;
    }
    public void setC(double[] c) {
        this.c = c;
    }
}