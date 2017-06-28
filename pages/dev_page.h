/*! \page page:dev For Developers
 * 
 * \tableofcontents
 * 
 * \section sec-interfaces Interfaces
 * 
 * \subsection MATLAB
 * 
 * To install <code>SuperSCS</code> in MATLAB, you need to build a MEX interface.
 * 
 * Do:
 * 
 * 
 * 
 * \code{.cpp}
 * cd matlab;
 * make_scs;
 * \endcode
 * 
 * 
 * \subsection CVX
 * 
 * Necessary steps:
 * 
 * - Download and unpack <code>CVX 3.0</code> from [here](http://cvxr.com/cvx/beta/)
 * - Navigate to the <code>cvx</code> directory
 * - Run <code>cvx_setup</code>
 * - To use <code>SuperSCS</code>, we call <code>cvx</code> with <code>cvx_solver scs</code>
 * and setting the parameter <code>do_super_scs</code> to 1.
 * 
 * Here is an example of an LP problem
 * ~~~~~{.m}
 * A(1,1) = 0.3; A(4,1) = -0.5;
 * A(2,2) = 0.7; A(4,2) =  0.9; A(3,3) = 0.2;
 * A = sparse(A);
 * n = size(A,2);
 * b = [0.2; 0.1; -0.1; 0.1];
 * c = [1;-2;-3];
 * 
 * cvx_begin
 *     cvx_solver scs
 *     cvx_solver_settings('eps', 1e-8, 'do_super_scs', 1, 'rho_x', 1,...
 *          'direction', 100, 'memory', 50);
 *     variable x(n);
 *     dual variable y;
 *     minimize( x'*x + c' * x );
 *     subject to
 *          y : A * x <= b;
 * cvx_end
 *~~~~~
 * 
 * We have chosen the SuperSCS mode with \f$\rho_x=1\f$, the restarted Broyden 
 * direction and memory equal to \f$50\f$.
 * 
 * We have set the tolerance to \f$10^{-8}\f$.
 * 
 * 
 */