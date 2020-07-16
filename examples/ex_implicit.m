%% Example of how to deal with implicitly define variables.

addpath('../');
warning('off','AutoDiff:autoext');

x = myAD(5);

y = sqrt(5);
y = implicit_pre(y, x);

eqns = y - sqrt(x);

y = implicit_post(y, eqns);
