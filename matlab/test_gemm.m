randn('seed',132);
m = 3;
n = 2;
k = 4;
A = randn(k,m);
B = randn(k,n);
C = randn(m,n);
alpha = randn;
beta = randn;

C2 = beta*C + alpha*A'*B;