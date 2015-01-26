clear all;close all;
randn('seed',1)
rand('seed',1)

N = 1e5;
test = 1000*randn(3,N);
a = rand(N,1); a(a>0.9) = 1; a(a<0.1) = 0;
for i=1:N
    [v(:,i), it(i)] = proj_pow(test(:,i),a(i));
    vd(:,i) = v(:,i) - test(:,i);
    err(i) = max((abs(v(3,i)) - (v(1,i)^a(i)) * (v(2,i)^(1-a(i))))/norm(test(:,i)), 0);
    errd(i) = max((abs(vd(3,i)) - ((vd(1,i)^a(i)) / (a(i)^a(i))) * ((vd(2,i)^(1-a(i)))/((1-a(i))^(1-a(i)))))/norm(test(:,i)), 0);
    ip(i) = abs(v(:,i)'*vd(:,i));
end
max(err)
max(errd)
max(ip)