%% Speed Test of Using Symbolic VS Automatic Differentiation
% by SeHyoun Ahn, Sept 2016
%%

%%
% An html file properly marked up comments is avaiable at <http://www.princeton.edu/~sehyouna/EXAMPLE_test_Syntax.html>

%%
% Evaluates the Jacobian of
%
% \$\$ y = \begin{bmatrix} x_1 & 0 & 0 & ... & 0 \\
%                     0 & x_2 & 0 & ... & 0\\
%                      0 & ... & ... &...& 0\\
%    0 & 0 & ... & 0 & x_n\end{bmatrix}
%     \begin{bmatrix} x_1\\ x_2\\ ...\\ x_n\end{bmatrix} \$\$
%     at \$(x_1, x_2 , ... , x_n) = (1,2,...,n)\$
%

n = 200;

%% Automatic Differentiation
fprintf('n = %d\n',n);
addpath('../');
fprintf('\n\nAUTOMATIC DIFFERENTIATION\n');
tAD = tic;
x = myAD((1:n)');
A = spdiags(x,0,n,n);
y = A*x;
dydx_val_AD = getderivs(y);
tAD = toc(tAD);
fprintf('   Done!\n\n');

%% Symbolic Differentiation
disp('SYMBOLIC')
tSymb = tic;

tSetup = tic;
vector = sym('x_%d',[1 n])';
A      = diag(vector);
y      = A*vector;
fprintf('  Initialize Symbolic Variables: %2.6f secs\n',toc(tSetup));

tDeriv = tic;
dydx   = jacobian(y,vector');
fprintf('  Time to Differentiate: %2.6f secs\n',toc(tDeriv));

tConvert  = tic;
dydx_eval = matlabFunction(dydx);
fprintf('  Time to Convert to Function: %2.6f secs\n',toc(tConvert));
fprintf('     You would only need to do this conversion once to evaluate \n');
fprintf('     derivatives at different values of x_n.\n');
tEval     = tic;
x = num2cell(1:n);
dydx_val = dydx_eval(x{:});
fprintf('  Time to Calculate Numerical Values: %2.6f secs\n\n\n',toc(tEval));
% Uncomment the following and comment above block to Test evaluation 
%     by using <subs>, but using matlabFunction is faster than subs
%{
tSubs    = tic;
dydx_val = subs(dydx,vector',1:n);
fprintf('Time to Calculate Numerical Values (Subs): %2.6f\n\n\n',toc(tSubs));
%}

tSymb = toc(tSymb);

fprintf('Total time using Automatic Differentiation: %2.6f secs\n',tAD);
fprintf('Total time using Symbolic: %2.6f secs\n',tSymb);
fprintf('   Note that just the evaluation of derivatives using symbolic\n');
fprintf('      is comparable to the entire computation of values and derivatives using AD.\n');
