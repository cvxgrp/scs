/*! \page page_benchmarks Benchmarks
 * 
 * \tableofcontents 
 * 
 * In what follows we give a few code snippets in MATLAB and compare SuperSCS
 * with SCS.
 * 
 * This is a preliminary study to show that SuperSCS outperforms SCS, but 
 * a more thorough analysis is necessary.
 * 
 * In what follows we compare the two algorithms using five different types of problems:
 * (i) a LASSO-like problem  (\f$\ell_1\f$-regularized least squares), (ii) 
 * a semidefinite program (SDP) and, in particular, a minimum-norm problem and an 
 * LMI-constrained problem, (iii) a logistic regression problem, (iv) a minimum 
 * \f$p\f$-norm problem, (v) a 2-norm-constrained minimum-norm problem and, last ,
 * (vi) a matrix completion problem involving the nuclear norm.
 * 
 * \note The following benchmarks are available in 
 * [/tests/matlab/](https://github.com/kul-forbes/scs/tree/master/tests/matlab).
 * 
 * \note Reported runtimes were recorded on an Intel® Core™ i5-6200U CPU @ 2.30GHz × 4
 * machine with 11.6 GiB RAM running 64-bit Ubuntu 14.04.
 * 
 * \note Both SuperSCS and SCS were compiled using <code>gcc</code> v4.8.5 and accessed
 * via MATLAB® using a MEX interface.
 * 
 * \remark For the sake of reproducibility, we use
 * 
 * ~~~~~{.m}
 * rng('default');
 * rng(1);
 * ~~~~~
 * 
 * \section sec_lasso LASSO-type problem
 * We solve a simple LASSO problem of the form
 * \f[
 * \mathrm{Minimize}_x\ \frac{1}{2} \|Ax-b\| + \mu \|x\|_1,
 * \f]
 * 
 * with \f$x\in\mathbb{R}^n\f$, \f$A\in\mathbb{R}^{m\times n}\f$.
 * 
 * Here \f$n=10000\f$, \f$m=2000\f$.
 * 
 * We define the problem data as follows:
 * 
 * ~~~~~{.m}
 * s = ceil(n/10);
 * x_true=[randn(s,1);zeros(n-s,1)];
 * x_true = x_true(randperm(n));
 * 
 * density = 0.1;
 * rcA = .1;
 * A=sprandn(m,n,density,rcA);
 * 
 * b = A*x_true + 0.1*randn(m,1);
 * mu = 1;
 * ~~~~~
 * 
 * Notice that \f$A\f$ is a sparse matrix with density \f$10\%\f$ and reciprocal
 * condition number \f$0.1\f$.
 * 
 * We then formulate the problem in CVX with accuracy \f$\epsilon=10^{-4}\f$.
 * 
 * We solve it using SuperSCS with \ref sec-restarted-broyden "restarted Broyden"
 * directions with memory \f$100\f$ and using \f$k_0=0\f$:
 * 
 * ~~~~~{.m}
 * do_super_scs = 0;
 * 
 * tic;
 * cvx_begin
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-4,...
 *         'do_super_scs',do_super_scs,...
 *         'direction', 100,...
 *         'memory', 100)
 *     variable x_c(n)
 *     minimize(0.5*sum_square(A*x_c - b) + mu*norm(x_c,1))
 * cvx_end
 * toc
 * ~~~~~
 * 
 * Notice that we may pass any parameters to SCS using CVX's <code>cvx_solver_settings</code>.
 * 
 * SuperSCS with the above options terminates after \f$169\f$ iterations and \f$85.1\f$ s.
 * 
 * Instead, SCS requires \f$3243\f$ iterations which corresponds to \f$802\f$ s.
 * 
 * SuperSCS exhibits a speed-up of \f$\times 9.4\f$.
 * 
 * \note Any values that are not overridden using <code>cvx_solver_settings</code> are 
 * set to their \ref ::setDefaultSettings "default values".
 * 
 * \section sec_SDP Semidefinite Programming
 * 
 * Here we solve the problem 
 * \f{eqnarray*}{
 *  &&\mathrm{Minimize}_{Z}\ \|Z-P\|_{\mathrm{fro}}\\
 *  &&Z\geq 0\\
 *  &&Z: \mathrm{Hermitian},\ \mathrm{Toeplitz}
 * \f}
 * ~~~~~{.m}
 * n=800;
 * P = randn(n,n);
 * cvx_begin sdp
 * cvx_solver scs
 * cvx_solver_settings('eps', 1e-4,...
 *         'do_super_scs',1,...
 *         'memory', 50)
 *     variable Z(n,n) hermitian toeplitz
 *     dual variable Q
 *     minimize( norm( Z - P, 'fro' ) )
 *     Z >= 0 : Q;
 * cvx_end
 * ~~~~~
 * 
 * SuperSCS solves this problem in \f$218\f$s (\f$103\f$ iterations), whereas
 * SCS requires \f$500\f$s to terminate (\f$392\f$ iterations).
 * 
 * 
 * Another interesting SDP problem is that of determining a symmetric positive 
 * definite matrix \f$P\f$ which solves
 * 
 * \f{eqnarray*}{
 *  &&\mathrm{Minimize}_{P}\ \mathrm{trace}(P)\\
 *  &&P=P'\\
 *  &&P \geq I\\
 *  &&A'P + PA \leq -I
 * \f}
 * 
 * For this example we chose a (stable) matrix \f$A\f$ with its eigenvalues 
 * uniformly logarithmically spaced between \f$-10^{-1}\f$ and \f$-10^{1}\f$.
 * 
 * ~~~~~{.m}
 * n=100;
 * A = diag(-logspace(-0.5, 1, n));
 * U = orth(randn(n,n));
 * A = U'*A*U;
 * 
 * tic
 * cvx_begin sdp
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-3,...
 *         'do_super_scs', 1, ...
 *         'memory', 100);
 *     variable P(n,n) symmetric
 *     minimize(trace(P))
 *     A'*P + P*A <= -eye(n)
 *     P >= eye(n)
 * cvx_end    
 * toc    
 * ~~~~~ 
 * 
 * SuperSCS converges after \f$40.7\f$s (\f$314\f$ iterations), whereas SCS converges
 * after \f$415\f$s (\f$4748\f$ iterations)
 * 
 * 
 * 
 * 
 * \section sec_log_reg Logistic Regression
 * 
 * Here we solve the following \f$\ell_1\f$-regularized logistic regression problem:
 * 
 * \f[
 * \mathrm{Minimize}_w\ \lambda \|w\|_1 - \sum_{i}\log(1+\exp(a' w_i + b))
 * \f]
 * 
 * ~~~~~{.m}
 * density = 0.1;
 * p = 1000;   % features
 * q = 10*p;   % total samples
 * w_true = sprandn(p,1,0.2);
 * X_tmp = 3*sprandn(p, q, density);
 * ips = -w_true'*X_tmp;
 * ps = (exp(ips)./(1 + exp(ips)))';
 * labels = 2*(rand(q,1) < ps) - 1;
 * X_pos = X_tmp(:,labels==1);
 * X_neg = X_tmp(:,labels==-1);
 * X = [X_pos -X_neg]; % include labels with data
 * lam = 2; 
 * ~~~~~
 * 
 * We formulate the problem as follows:
 * 
 * ~~~~~{.m}
 * tolerance = 1e-3;
 * cvx_begin
 *    cvx_solver scs
 *    cvx_solver_settings('eps', tolerance,...
 *      'do_super_scs', 1,...
 *      'ls', 5,...
 *      'memory', 50)
 *    variable w(p)
 *    minimize(sum(log_sum_exp([sparse(1,q);w'*X])) + lam * norm(w,1))       
 * cvx_end
 * ~~~~~
 * 
 * For \f$\epsilon=10^{-3}\f$, SuperSCS terminates after \f$131\f$s (\f$134\f$ iterations), 
 * while SCS requires \f$329\f$s (\f$675\f$ iterations).
 * 
 * 
 * For \f$\epsilon=10^{-4}\f$, SuperSCS requires \f$179\f$s (\f$183\f$ iterations),
 * whereas SCS still requires as much as \f$497\f$s (\f$983\f$ iterations).
 * 
 * In both cases, SuperSCS is about \f$2.5\f$ times faster.
 * 
 * 
 * 
 * 
 * \section sec_p_norm Minimization of p-norm
 * 
 * Consider the problem
 * \f{eqnarray*}{
 *   &&\mathrm{Minimize}_x\ \|x\|_p\\
 *   &&Gx = f
 * \f}
 * 
 * This is the problem formualation in CVX
 * 
 * ~~~~~{.m}
 * n = 2000;
 * m = ceil(n/4);
 * 
 * density = 0.1;
 * G = sprandn(m,n,density);
 * f = randn(m,1) * n * density;
 * pow = 1.5;
 * 
 * tic
 * cvx_begin
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-4,...
 *         'do_super_scs', 0, ...
 *         'memory', 50,...
 *         'use_indirect', 0);
 *     variable x(n)
 *     minimize(norm(x, pow))
 *     subject to
 *       G*x == f
 *     cvx_end
 * toc
 * ~~~~~
 * 
 * For \f$\epsilon=10^{-3}\f$, SuperSCS converged after \f$11.3\f$s (\f$1355\f$ iterations) ,
 * while SCS required \f$39\f$s (\f$4911\f$ iters).
 * 
 * For \f$\epsilon=10^{-4}\f$, SuperSCS terminated after 17s (\f$4909\f$ iterations),
 * while SCS failed to terminate within \f$10^4\f$ iterations.
 * 
 * 
 * 
 * \section sec_another_problem Norm-constrained minimum-norm problem
 * 
 * Here we solve a constrained problem of the form
 * 
 * \f{eqnarray*}{
 * &&\mathrm{Minimize}_x\ \|Ax-b\|\\
 * &&\|Gx\| \leq 1
 * \f}
 * 
 * with \f$x\in\mathbb{R}^m\f$, \f$A\in\mathbb{R}^{m\times n}\f$ and 
 * \f$G\in\mathbb{R}^{2n\times n}\f$.
 * 
 * The problem data are:
 * 
 * ~~~~~{.m}
 * m = 3000; 
 * n = ceil(m/2);
 * A = sprandn(m,n,0.5,0.01);
 * b = 10*randn(m,1);
 * G = 2*sprandn(2*n, n, 0.1);
 * ~~~~~
 * 
 * And we need to solve it up to an accuracy of \f$\epsilon=10^{-3}\f$.
 * 
 * The problem is formulated in CVX as follows:
 * 
 * ~~~~~{.m}
 * tic;
 * cvx_begin
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-3,...
 *         'do_super_scs', 1, ...
 *         'memory', 20);
 *     variable x(n)
 *     minimize( norm(A*x-b) )
 *     subject to
 *         norm(G*x) <= 1
 * cvx_end
 * toc
 * ~~~~~
 * 
 * SuperSCS converges after \f$37.1\f$s (\f$116\f$ iterations).
 * 
 * On the other hand, SCS did not converge after \f$242s\f$ (\f$5000\f$ iterations).
 * 
 * Note that SuperSCS can attain much higher precision; for instance, 
 * it converges with \f$\epsilon=10^{-4}\f$ in \f$41.2\f$s (\f$131\f$ iterations) and
 * with \f$\epsilon=10^{-6}\f$ it converges after \f$58\f$s (\f$194\f$ iterations). 
 * 
 * 
 * \section sec_matrix_completion Regularized matrix completion problem
 * 
 * Let \f$M\f$ be a given matrix whose elements \f$\{(i,j)\}_{i\in I, j\in J}\f$
 * are missing.
 * 
 * Here, we formualte the matrix completion problem as a nuclear norm minimization 
 * problem.
 * 
 * \f{eqnarray*}{
 * &&\mathrm{Minimize}_{X}\ \|X-M\|_* + \lambda \|X\|_{\mathrm{fro}}^2\\
 * &&X_{i',j'} = M_{i',j'},\ \forall i'\notin I,\ j'\notin J
 * \f}
 * 
 * 
 * 
 * ~~~~~{.m}
 * rng('default');
 * rng(1);
 * m=200;
 * n=200;
 * n_nan = ceil(0.8*m*n);
 * 
 * M = sprandn(m, n, 0.4);
 * idx = randperm(m*n);
 * M(idx(1:n_nan))=nan;
 * lam = 0.5;
 * 
 * tic
 * cvx_begin sdp
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-3,...
 *         'do_super_scs', 1,...
 *         'memory', 100)
 *     variable X(m,n)
 *     minimize (norm_nuc(X)  + lam*sum_square(X(:)))
 *     subject to
 *     for i=1:m
 *         for j=1:n
 *             if (~isnan(M(i,j)))
 *                 X(i,j)==M(i,j)
 *             end
 *         end
 *     end
 * cvx_end
 * toc
 * ~~~~~
 * 
 * SuperSCS converges in \f$18.4s\f$ (\f$102\f$ iterations), whereas SCS takes
 * \f$270s\f$ (\f$6061\f$ iterations).
 * 
 * Changing <code>lam</code> to <code>0.05</code> deems SCS unable to converge 
 * within 10,000 iterations (gain, with \f$\epsilon=10^{-3}\f$), while SuperSCS converges 
 * in \f$26.2\f$s (\f$145\f$ iterations).
 * 
 * 
 * \section sec_portfolio Portfolio selection
 * 
 * Here we solve the following portfolio selection problem:
 * 
 * \f{eqnarray*}{
 * &&\mathrm{Maximize}_z\ \mu'z  - \gamma z'\Sigma z\\
 * &&1'z = 1\\
 * &&z \geq 0,
 * \f}
 * 
 * where \f$z\f$ is the portfolio of assets, \f$\mu\f$ is the vector of 
 * expected returns, \f$\gamma\f$ is a risk-aversion parameter and 
 * \f$\Sigma=F'F+D\f$ is the asset return covariance matrix with \f$F\f$
 * being the <em>factor loading matrix</em> and \f$D\f$ is a diagonal matrix 
 * called the <em>asset-specific risk</em>.
 * 
 * The problem is generated as follows: 
 * 
 * ~~~~~{.m}
 * density = 0.1;
 * rc = 0.5; % estimated reciprocal condition number
 * 
 * n = 100000;
 * m = 100;
 * 
 * mu = exp(0.01 * randn(n, 1)) - 1; % returns
 * D = rand(n,1) / 10; % idiosyncratic risk
 * F = sprandn(n, m, density, rc) / 10; % factor model
 * gamma = 1;
 * B = 1;
 * ~~~~~
 * 
 * The problem is solved using CVX as follows:
 * 
 * ~~~~~{.m}
 * cvx_begin
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-4,...
 *         'do_super_scs', 1,...
 *         'direction', 100,...
 *         'memory', 100)
 *     variable x(n)
 *     maximize(mu'*x - gamma*(sum_square(F'*x) + sum_square(D.*x)))r
 *     sum(x) == B
 *     x >= 0
 * cvx_end
 * ~~~~~
 * 
 * The above problem is solved in \f$46.7\f$s and at \f$292\f$ iterations with SuperSCS 
 * running with restarted Broyden directions with memory equal to \f$100\f$.
 * 
 * The respective results for SCS are \f$157\f$s and \f$3050\f$ iterations (approximately
 * \f$3.3\f$ times slower).
 * 
 * For a lower tolerance of \f$\epsilon=10^{-3}\f$ , SuperSCS with memory equal to \f$100\f$
 * terminates at \f$162\f$ iterations after \f$23.5\f$s while SCS converges after \f$43.9\f$s 
 * at \f$787\f$ iterations.
 * 
 * 
 * \section sec_pca Principal components analysis 
 * 
 * Here we solve the following PCA problem (using an \f$\ell_1\f$-regularization).
 * 
 * The problem has the form
 * 
 * \f{eqnarray*}{
 *  &&\mathrm{Maximize}_X\ \mathrm{trace}(SX) - \lambda \|X\|_1\\
 * && \mathrm{trace}(X) = 1\\
 * && X = X' \succeq 0
 * \f}
 * 
 * ~~~~~{.m}
 * d = 200;
 * p = 10;
 * A = sprandn(p,d,0.3,0.0001);
 * S = full(A'*A);
 * lam = 2;
 * 
 * cvx_begin sdp
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-3,...
 *        'verbose', 1,...
 *        'do_super_scs', 0, ...
 *        'direction', 100, ...
 *        'memory', 100);
 *     variable X(d,d) symmetric
 *     minimize(-trace(S*X) + lam*norm(X,1))
 *     trace(X)==1
 *     X>=0
 * cvx_end  
 * ~~~~~
 * 
 * SuperSCS: \f$7.6\f$s, \f$95\f$ iterations
 * 
 * SCS: \f$65\f$s, \f$2291\f$ iterations
 * 
 * Speed-up: \f$\times 8.6\f$
 */