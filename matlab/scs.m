function [ varargout ] = scs( varargin )
% scs 1.0
data = varargin{1};
K = varargin{2};
pars = varargin{3};
if (isfield(pars,'USE_INDIRECT') && pars.USE_INDIRECT)
    [  varargout{1:nargout}  ] = scs_indirect( data, K, pars);
else
    [  varargout{1:nargout}  ] = scs_direct( data, K, pars);
end