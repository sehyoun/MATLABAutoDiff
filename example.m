%% Examples of using automatic differentiation
% Written by SeHyoun Ahn, June 2016

%%
% f_1(x,y,z) = x^2+cos(y)*(x+z)+exp(z) at (x,y,z)=(1,2,3)
% f_2(x,y,z) = yz;

v = myAD([1;2;3]);
x=v(1); y=v(2); z=v(3);
f=[x^2+cos(y)*(x+z)+exp(z);
    y*z];
disp('The function evaluated at (1,2,3) is');
disp(getvalues(f));
disp('the derivatives evaluated at (1,2,3) are');
disp(full(getderivs(f)));

%%
% An example of matrix vector multiplication
%
%     [ x y 0 ]   [ x ]   [ x^2+y^2    ]
% f = [ y y z ] * [ y ] = [ xy+y^2+z^2 ]
%     [ 0 z z ]   [ z ]   [ yz+z^2     ]
% evaluated at (x,y,z)=(1,2,3)

v=myAD([1:3]');
A=spdiags(v,0,3,3)+spdiags(v(2:3),-1,3,3)+spdiags(v,1,3,3);
f=A*v;
disp('f is');
disp(getvalues(f));
disp('The Jacobian of f evaluated at (1,2,3) is')
disp(full(getderivs(f)));

%%
% An example of using fsolve.
% x^2+y^4+x*(z-0.2)*e^z=5 solving for dz/dx and dz/dy at (x,y)=(1,1)

% Refer to README.pdf for syntactic requirement on definition of functions.
% x = v(2), y = v(3), z = v(1)
func = @(v) v(2).^2+v(3).^4+v(2)*(v(1)-0.2)*exp(v(1))-5;
v0 = myAD([1;1]);            % Initialize x and y
z=fsolve(func, 0.5 ,v0);     % Refer to documentation for syntax
disp('Solution of z is given by');
disp(getvalues(z));
disp('The derivatives at (x,y)=(1,1) are');
disp(getderivs(z));
