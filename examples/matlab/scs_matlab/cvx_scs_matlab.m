function shim = cvx_scs_matlab( shim )

% CVX_SOLVER_SHIM	SCS interface for CVX.
%   This procedure returns a 'shim': a structure containing the necessary
%   information CVX needs to use this solver in its modeling framework.

if ~isempty( shim.solve ),
    return
end
if isempty( shim.name ),
    fname = 'scs_matlab.m';
    flen = length(fname);
    shim.name = 'SCS_MATLAB';
    shim.config = struct( 'dualize', 1, 'nonnegative', 1, 'lorentz', 1, ...
        'semidefinite', 1, 'exponential', 1 );
    shim.warning = 'The algorithm employed by SCS typically requires a large number of iterations to obtain the levels of accuracy normally sought by CVX. Consider reducing the precision using the CVX_PRECISION command to improve performance.';
    shim.params.newsdp = false;
    fpaths = which( fname, '-all' );    
    if ~iscell(fpaths),
        fpaths = { fpaths };
    end
    old_dir = pwd;
    oshim = shim;
    shim = [];
    for k = 1 : length(fpaths),
        fpath = fpaths{k};
        if ~exist( fpath, 'file' ) || any( strcmp( fpath, fpaths(1:k-1) ) ),
            continue
        end
        new_dir = fpath(1:end-flen-1);
        cd( new_dir );
        tshim = oshim;
        tshim.fullpath = fpath;
        tshim.version = 'unknown';
        tshim.location = new_dir;
        otp = [];
        try
            if exist( [ fpath(1:end-2), '_version.', mexext ], 'file' ),
                otp = scs_version;
            end
        catch errmsg
            tshim.error = sprintf( 'Unexpected error:\n%s\n', errmsg.message );
        end
        if isempty(otp) && isempty(tshim.error),
            try
                fid = fopen(fname);
                otp = fread(fid,Inf,'uint8=>char')';
                fclose(fid);
                otp = regexp( otp, 'scs \d+\.\d+\.?\d+?', 'match' );
                if ~isempty(otp), 
                    otp = otp{1}(5:end);
                end
            catch errmsg
                tshim.error = sprintf( 'Unexpected error:\n%s\n', errmsg.message );
            end
        end
        if isempty(tshim.error),
            if ~isempty(otp), 
                tshim.version = otp; 
                otp = otp(1:find(otp=='.',1,'last')-1);
                tshim.eargs.newsdp = str2double(otp) >= 1.1;
            else
                tshim.eargs.newsdp = false;
            end
            tshim.params.newsdp = tshim.eargs.newsdp;
            tshim.solve = @solve;
            if k ~= 1,
                tshim.path = [ new_dir, pathsep ];
            end
        end
        shim = [ shim, tshim ]; %#ok
    end
    cd( old_dir );
    if isempty( shim ),
        shim = oshim;
        shim.error = 'https://github.com/cvxgrp/scs';
    end
else
    shim.solve = @solve;
    if isfield( shim.params, 'newsdp' ),
        shim.eargs.newsdp = shim.params.newsdp;
    else
        shim.eargs.newsdp = false;
    end
end

function [ x, status, tol, iters, y ] = solve( At, b, c, nonls, params )

n = length( c );
m = length( b );
K = struct( 'f', 0, 'l', 0, 'e', 0, 'q', [], 's', [] );
reord = struct( 'n', 0, 'r', [], 'c', [], 'v', [] );
reord = struct( 'f', reord, 'l', reord, 'e', reord, 'q', reord, 's', reord );
reord.f.n = n;
for k = 1 : length( nonls ),
    temp = nonls( k ).indices;
    nn = size( temp, 1 );
    nv = size( temp, 2 );
    nnv = nn * nv;
    tt = nonls( k ).type;
    reord.f.n = reord.f.n - nnv;
    switch tt,
        case 'nonnegative',
            reord.l.r = [ reord.l.r ; temp(:) ];
            reord.l.c = [ reord.l.c ; reord.l.n + ( 1 : nnv )' ];
            reord.l.v = [ reord.l.v ; ones( nnv, 1 ) ];
            reord.l.n = reord.l.n + nnv;
        case 'exponential',
            reord.e.r = [ reord.e.r ; temp(:) ];
            reord.e.c = [ reord.e.c ; reord.e.n + ( 1 : nnv )' ];
            reord.e.v = [ reord.e.v ; ones(nnv,1) ];
            reord.e.n = reord.e.n + nnv;
            K.e = K.e + nv;
        case 'lorentz',
            temp = temp( [ end, 1 : end - 1 ], : );
            reord.q.r = [ reord.q.r ; temp(:) ];
            reord.q.c = [ reord.q.c ; reord.q.n + ( 1 : nnv )' ];
            reord.q.v = [ reord.q.v ; ones(nnv,1) ];
            reord.q.n = reord.q.n + nnv;
            K.q = [ K.q, nn * ones( 1, nv ) ];
        case 'semidefinite',
            n2 = 0.5 * ( sqrt( 8 * nn + 1 ) - 1 );
            K.s = [ K.s, n2 * ones( 1, nv ) ];
            if true, %params.newsdp,
                scale = ones(nn,nv);
                sdiag = cumsum([1,n2:-1:2]);
                scale(sdiag,:) = sqrt(2);
                reord.s.r = [ reord.s.r ; temp(:) ];
                reord.s.c = [ reord.s.c ; reord.s.n + ( 1 : nnv )' ];
                reord.s.v = [ reord.s.v ; scale(:) ];
                reord.s.n = reord.s.n + nnv;
            else
                str = cvx_create_structure( [ n2, n2, nv ], 'symmetric' );
                [ rr, cc, vv ] = find( cvx_invert_structure( str, true ) );
                rr = temp( rr );
                reord.s.r = [ reord.s.r; rr( : ) ];
                reord.s.c = [ reord.s.c; cc( : ) + reord.s.n ];
                reord.s.v = [ reord.s.v; vv( : ) ];
                reord.s.n = reord.s.n + n2 * n2 * nv;
            end
    end
end
if reord.f.n > 0,
    reord.f.r = ( 1 : n )';
    reord.f.r( [ reord.l.r ; reord.e.r ; reord.q.r ; reord.s.r ] ) = [];
    reord.f.c = ( 1 : reord.f.n )';
    reord.f.v = ones(reord.f.n,1);
end
n_d = max( m - n - reord.f.n + 1, isempty( At ) );
if n_d,
    reord.l.n = reord.l.n + n_d;
end
K.f = reord.f.n;
K.l = reord.l.n;
                               n_out =         reord.f.n;
reord.l.c = reord.l.c + n_out; n_out = n_out + reord.l.n;
reord.q.c = reord.q.c + n_out; n_out = n_out + reord.q.n;
reord.s.c = reord.s.c + n_out; n_out = n_out + reord.s.n;
reord.e.c = reord.e.c + n_out; n_out = n_out + reord.e.n;
reord = sparse( ...
    [ reord.f.r ; reord.l.r ; reord.q.r ; reord.s.r ; reord.e.r ], ...
    [ reord.f.c ; reord.l.c ; reord.q.c ; reord.s.c ; reord.e.c ], ...
    [ reord.f.v ; reord.l.v ; reord.q.v ; reord.s.v ; reord.e.v ], ...
    n, n_out );

At = reord' * At;
c  = reord' * c;

data.A = sparse(At);
data.b = full(c);
data.c = -full(b);

KK = K;
K  = struct;
if KK.f, K.f = KK.f; end
if KK.l, K.l = KK.l; end
if ~isempty( KK.q ), K.q = KK.q(:); end
if ~isempty( KK.s ), K.s = KK.s(:); end
if ~isempty( KK.e ), K.ed = KK.e(:); end
prec = params.precision;
pars.verbose = +(~params.quiet);
pars.max_iters = 10000;
pars.eps = max(eps,prec(1));

[ yy, xx, ~, info ] = cvx_run_solver( @scs_matlab, data, K, pars, 'xx', 'yy', 'info', 3, params );

xx = full( xx );
yy = full( yy );
x  = real( reord * xx );
y  = yy;

iters = info.iter;
switch info.status
    case { 'Solved', 'Solved/Inaccurate', 'Inaccurate/Solved' },
        tol = max([info.resPri,info.resDual,info.relGap]);
        status = 'Solved';
    case { 'Unbounded', 'Unbounded/Inaccurate', 'Inaccurate/Unbounded' },
        tol = info.resPri;
        status = 'Infeasible';
    case { 'Infeasible', 'Infeasible/Inaccurate', 'Inaccurate/Infeasible' },
        status = 'Unbounded';
        tol = info.resDual;
    otherwise,
        tol = Inf;
        status = 'Failed';
end
if tol > prec(2),
    if tol > prec(3)
        status = 'Failed';
    elseif status(1) ~= 'F',
        status = [ 'Inaccurate/', status ];
    end
end

% Copyright 2012 CVX Research, Inc.
% See the file COPYING.txt for full copyright information.
% The command 'cvx_where' will show where this file is located.
