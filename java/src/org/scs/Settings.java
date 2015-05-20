package org.scs;

public class Settings {
    // set to defaults:
    private int maxIters = 2500;        /* maximum iterations to take: 2500 */
    private double eps = 1e-3;          /* convergence tolerance: 1e-3 */
    private double alpha = 1.8;         /* relaxation parameter: 1.8 */
    private double rhoX = 1e-3;;        /* x equality constraint scaling: 1e-3 */
    private double cgRate = 2;          /* for indirect, tolerance goes down like (1/iter)^cg_rate: 2 */
    private boolean verbose = true;     /* boolean, write out progress: 1 */
    private boolean normalize = true;   /* boolean, heuristic data rescaling: 1 */
    private double scale = 5;           /* if normalized, rescales by this factor: 5 */
    private boolean warmStart = false;  /* boolean, warm start (put initial guess in Sol struct): 0 */

    // getters:
    public int getMaxIters() {
        return maxIters;
    }
    public double getEps() {
        return eps;
    }
    public double getAlpha() {
        return alpha;
    }
    public double getRhoX() {
        return rhoX;
    }
    public double getCgRate() {
        return cgRate;
    }
    public boolean isVerbose() {
        return verbose;
    }
    public boolean isNormalize() {
        return normalize;
    }
    public double getScale() {
        return scale;
    }
    public boolean isWarmStart() {
        return warmStart;
    }

    // setters:
    public void setMaxIters(int maxIters) {
        this.maxIters = maxIters;
    }
    public void setEps(double eps) {
        this.eps = eps;
    }
    public void setAlpha(double alpha) {
        this.alpha = alpha;
    }
    public void setRhoX(double rhoX) {
        this.rhoX = rhoX;
    }
    public void setCgRate(double cgRate) {
        this.cgRate = cgRate;
    }
    public void setVerbose(boolean verbose) {
        this.verbose = verbose;
    }
    public void setNormalize(boolean normalize) {
        this.normalize = normalize;
    }
    public void setScale(double scale) {
        this.scale = scale;
    }
    public void setWarmStart(boolean warmStart) {
        this.warmStart = warmStart;
    }
}
