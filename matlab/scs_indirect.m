function [x, y, s, info] = scs_indirect(data, cone, params)			    
% Operator-splitting method for solving cone problems (indirect)
%
% This implements a cone solver. It solves:
%
% min. c'x
% subject to Ax + s = b
% s \in K
%
% where x \in R^n, s \in R^m
%
% this uses the indirect linear equation solver version of SCS
%
% K is product of cones in this particular order:
% free cone, lp cone, second order cone(s), semi-definite cone(s), primal
% exponential cones, dual exponential cones
%
% data must consist of data.A, data.b, data.c, where A,b,c used as above.
%  
% cone struct must consist of:
% cone.f, length of free cone (for equality constraints)
% cone.l, length of lp cone
% cone.q, array of SOC lengths
% cone.s, array of SD lengths
% cone.ep, number of primal exp cones
% cone.ed, number of dual exp cones
%
% Optional fields in the params struct are:
%   alpha       : over-relaxation parameter, between (0,2).
%   rho_x       : momentum of x term (1e-3 works well)
%   max_iters   : maximum number of ADMM iterations.
%   eps         : accuracy of solution
%   verbose     : verbosity level (0 or 1)
%   normalize   : heuristic data rescaling (0 or 1, off or on)
%   scale       : rescales data up by this factor (only used if normalize=1)
%   cg_rate     : the rate at which the CG tolerance is tightened (higher is tighter)
%
% to warm-start the solver add guesses for (x, y, s) to the data struct
%
error ('scs_indirect mexFunction not found') ;
