%% Test the differentiation of the SVD
%
% by SeHyoun Ahn, July 2019
%

n_div = 5;
dt = zeros(n_div,1);
dt(1) = 1e-5;

x = myAD(randn(n_div,1));

A_func = @(x) [x(1), x(2), x(4); 0, x(3), x(5)];
A = A_func(x);

[u, s, v] = svd(full(getvalues(A)), 'econ');

[u_diff, s_diff, v_diff] = svd_helper(A, u, s, v);

u_dt = u(:) + getderivs(u_diff)*dt;
s_dt = s(:) + getderivs(s_diff)*dt;
v_dt = v(:) + getderivs(v_diff)*dt;

x = x + dt;
A = A_func(x);
[u_new, s_new, v_new] = svd(full(getvalues(A)), 'econ');

max(abs(u_new(:) - u_dt))
max(abs(s_new(:) - s_dt))
max(abs(v_new(:) - v_dt))
