function [ x, y, info ] = scs( varargin )
% scs 1.0
data = varargin{1};
K = varargin{2};
if nargin >= 3
    pars = varargin{3};
else
    pars = [];
end
if (isfield(pars,'USE_INDIRECT') && pars.USE_INDIRECT)
    [  x, y, s, info  ] = scs_indirect( data, K, pars);
else
    [  x, y, s, info  ] = scs_direct( data, K, pars);
end