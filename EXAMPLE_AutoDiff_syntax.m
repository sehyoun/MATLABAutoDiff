%% Examples of Using Automatic Differentiation
% by SeHyoun Ahn, June 2016
%%

%%
% An html file properly marked up comments is avaiable at <http://www.princeton.edu/~sehyouna/EXAMPLE_AutoDiff_Syntax.html>

% Add folder containing <@myAD>. In this case, it is in the parent folder.
addpath('../');

%% Simple Example
% \$\$f_1(x,y,z) = x^2+cos(y)*(x+z)+e^z\$\$
% \$\$f_2(x,y,z) = yz\$\$
% at \$(x,y,z) = (1,2,3)\$

v = myAD([1;2;3]);
x=v(1); y=v(2); z=v(3);
f=[x^2+cos(y)*(x+z)+exp(z);
    y*z];
disp('The function evaluated at (1,2,3) is');
disp(getvalues(f));
disp('the derivatives evaluated at (1,2,3) are');
disp(full(getderivs(f)));

%% Example of Matrix Vector Multiplication
% \$\$ f = \begin{bmatrix} x & y & 0 \\ y& y& z \\ 0 & z & z
% \end{bmatrix} \begin{bmatrix}x\\y\\z \end{bmatrix} = 
% \begin{bmatrix} x^2+y^2\\  xy+y^2+z^2 \\ yz+z^2\end{bmatrix} \$\$
% evaluated at \$(x,y,z)=(1,2,3)\$

v=myAD([1:3]');
A=spdiags(v,0,3,3)+spdiags(v(2:3),-1,3,3)+spdiags(v,1,3,3);
f=A*v;
disp('f is');
disp(getvalues(f));
disp('The Jacobian of f evaluated at (1,2,3) is')
disp(full(getderivs(f)));

%%
% If you plan to use matrix multiplication in high dimensions,
% refer to the documentation to compile c-files for speed gains


%% Example of Using fsolve
% Given an implicitly defined variable \$z\$ from the relationship
% \$\$x^2+y^4+x\cdot (z-0.2)\cdot e^z=5 \$\$ solving for \$\frac{dz}{dx}\$ and \$\frac{dz}{dy}\$ evaluated at \$(x,y)=(1,1)\$

% Refer to README.pdf for syntactic requirement on definition of functions.
% x = v(2), y = v(3), z = v(1)
func = @(v) v(2).^2+v(3).^4+v(2)*(v(1)-0.2)*exp(v(1))-5;
v0 = myAD([1;1]);            % Initialize x and y
z=fsolve(func, 0.5 ,v0);     % Refer to documentation for syntax
disp('Solution of z is given by');
disp(getvalues(z));
disp('The derivatives at (x,y)=(1,1) are');
disp(getderivs(z));
