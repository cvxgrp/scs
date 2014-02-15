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
% cone struct is only used in cones.c, to add other cones
% simply add the relevant size data to the cone struct and edit the
% cones.c to include projection onto the new cone AND
% the dual cone
%
% Optional fields in the params struct are:
%   ALPHA       : over-relaxation parameter, between (0,2).
%   RHO_X       : momentum of x term (1e-3 works well)
%   MAX_ITERS   : maximum number of ADMM iterations.
%   EPS     : accuracy of solution
%   VERBOSE     : verbosity level (0 or 1)
%   NORMALIZE   : heuristic data rescaling (0 or 1, off or on)
error ('scs_indirect mexFunction not found') ;
