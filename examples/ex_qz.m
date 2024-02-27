%% Test the differentiation of the QZ decomposition
%
% by SeHyoun Ahn, July 2019
%
% Warning:
%   QZ decomposition is not unique up to rotations, so make sure that it is
%   testing two decompsitions with same ordering. This is the most likely case
%   with the <QZ> method of MATLAB.
%

n_div = 3;
dt = zeros(n_div, 1);
dt(1) = 1e-4;

addpath('../');
x = myAD(randn(n_div,1));

A_func = @(x) [x(1), x(2); x(3), x(3)];
B_func = @(x) [x(2), x(3); x(1), x(2)];

A = A_func(x);
B = B_func(x);

[AA, BB, Q, Z] = qz(full(getvalues(A)), full(getvalues(B)), 'complex');
[AA_diff, BB_diff, Q_diff, Z_diff] = qz_helper(A, B, AA, BB, Q, Z);

AA_dt = AA(:) + getderivs(AA_diff)*dt;
BB_dt = BB(:) + getderivs(BB_diff)*dt;
Q_dt = Q(:) + getderivs(Q_diff)*dt;
Z_dt = Z(:) + getderivs(Z_diff)*dt;

x = x + dt;
A = A_func(x);
B = B_func(x);

[AA_new, BB_new, Q_new, Z_new] = qz(full(getvalues(A)), full(getvalues(B)), 'complex');

max(abs(AA_new(:) - AA_dt))
max(abs(BB_new(:) - BB_dt))
max(abs(Q_new(:) - Q_dt))
max(abs(Z_new(:) - Z_dt))
