package org.scs;

public class Cone {
    private int f;
    private int l;
    private int[] q;
    private int[] s;
    private int ep;
    private int ed;
    private double[] p;

    // getters:
    public int getF() {
        return f;
    }
    public int getL() {
        return l;
    }
    public int[] getQ() {
        return q;
    }
    public int[] getS() {
        return s;
    }
    public int getEp() {
        return ep;
    }
    public int getEd() {
        return ed;
    }
    public double[] getP() {
        return p;
    }

    // setters:
    public void setF(int f) {
        this.f = f;
    }
    public void setL(int l) {
        this.l = l;
    }
    public void setQ(int[] q) {
        this.q = q;
    }
    public void setS(int[] s) {
        this.s = s;
    }
    public void setEp(int ep) {
        this.ep = ep;
    }
    public void setEd(int ed) {
        this.ed = ed;
    }
    public void setP(double[] p) {
        this.p = p;
    }

}