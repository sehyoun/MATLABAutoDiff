%% Test the differentiation of the Schur decomposition
%
% by SeHyoun Ahn, July 2019
%
% Warning:
%   Schur decomposition is not unique up to rotations, so make sure that it is
%   testing two decompsitions with same ordering. This is the most likely case
%   with the <schur> method of MATLAB. Also, with complex decompositions,
%   the eigenvector can be negative.
%

n_div = 25;
dt = zeros(n_div, 1);
dt(2) = 0.0001;

addpath('../');
x = myAD(randn(n_div,1));

A_func = @(x) reshape(x,5,5);
A = A_func(x);

[u, t] = schur(full(getvalues(A)), 'complex');
[u, t] = ordschur(u, t, 'lhp');
[u_diff, t_diff] = schur_helper(A, u, t);

u_dt = u(:) + getderivs(u_diff)*dt;
t_dt = t(:) + getderivs(t_diff)*dt;

x = x + dt;
A = A_func(x);
[u_new, t_new] = schur(full(getvalues(A)), 'complex');
[u_new, t_new] = ordschur(u_new, t_new, 'lhp');

max(abs(u_new(:) - u_dt))
max(abs(t_new(:) - t_dt))
