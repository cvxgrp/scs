clear;
rng(0)

l = 3;
tau = 0.21;
S=[];
U=[];
sig=[];
iter=1;

optsDir.memory = 5;
optsDir.delta = 0.5;
for i=1:15,
    Rx = randn(l,1);
    Rxold = randn(l,1);
    Sk = randn(l,1);
    Yk = randn(l,1);
    [d,S,U] = ls_broyden(Rx,Rxold,tau,Sk,Yk,S,U,sig,i,optsDir);    
end