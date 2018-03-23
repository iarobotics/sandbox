%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% matrix_hw_gr671 %
% Authors: Amalia-Lelia Cretu (ROB6)
%          Ion Sircu (ROB6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Exercise 6.1
close all
clc
clear all

cvx_begin
    variables x y
    minimise((x^2)+y+4)
    x-y-6>=0
    -x^2-(y+4)^2+16>=0
cvx_end

%[X,Y] = meshgrid(-10:.1:10,-10:.1:10);
x = -10:.2:10;
y = -10:.2:10;
[X,Y] = meshgrid(x,y)

F = X.^2+Y+4;
C2 = X-Y-6;


%Z3 = X.^3+Y.^3;

figure(1)
[C,h] = contour(X,Y,F);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)
hold on
[C,h] = contour(X,Y,C2);
set(h,'ShowText','on','TextStep',get(h,'LevelStep')*2)

%surf(X,Y,F)
%hold on
%surf(X,Y,C1)

hold on

%% Exercise 6.3
close all
clc
clear all

% Calculating gradient and Hessian in syms
syms u t
F = 1/4*((u^2)+4*(t^2)-4*(3*u+8*t)+100);
g=gradient(F,[u, t]);
H=hessian(F,[u, t]);

% Solving for the minimizer using CVX
cvx_begin quiet
    variables x y
    Z = 1/4*((x^2)+4*(y^2)-4*(3*x+8*y)+100);
    minimize(Z)
    x==2
    y>=0
cvx_end

% Calculating gradient and Hessian at minimizer
g_xy=double(subs(g, {u,t}, {x,y}))
H_xy=double(subs(H, {u,t}, {x,y}))

% Verifying second order derivative conditions
syms d d2
d=[0; d2];

%first_condition=g_xy*d




%% Exercise 8.1.c
close all
clc
clear all

A = [1 -1; 3 3];
b = [-2; -3];

cvx_begin
   variable x(2)
   variable l(2,1) % l for lambda
   minimize(x(1)^2 + x(2)^2 +5)
   subject to
     2*x - A'*l == 0;
     x(1)-x(2)+2 == 0;
     3*x(1)+3*x(2)+3 == 0;
cvx_end

%% Exercise 8.1.d
close all
clc
clear all

C = [1 -1; 3 3];
d = [-2; -3];

cvx_begin
   variable x(2)
   variable m(2,1) % m for miu
   minimize(x(1)^2 + x(2)^2 +5)
   subject to
     2*x - C'*m == 0;
     C'*m == 0;
     C >= 0;
     m >= 0;
     x(1)-x(2)+2 >= 0;
     3*x(1)+3*x(2)+3 >= 0;
cvx_end

%% Exercise 8.1.e
close all
clc
clear all

C = [3 3];
d = -3;
A = [1 -1];
b = -2;

cvx_begin
   variables x(2) l(1,1) m(1,1) % l for lambda, m for miu
   minimize(x(1)^2 + x(2)^2 +5)
   subject to
     2*x - A'*l - C'*m == 0;
     A*x == b;
     C'*m == 0;
     C >= 0;
     m >= 0;
     x(1)-x(2)+2 >= 0;
     3*x(1)+3*x(2)+3 >= 0;
cvx_end

%% Exercise 9.1.a
close all
clc
clear all

CT = [-1 -2 3];
C=transpose(CT);
A = [1 1 1];
b = 1;
D = [1 0 0; 0 1 0; 0 0 1]

for t=1:-.1:0
cvx_begin sdp
   variables x(3) L(1,1) m(3,1) % l for lambda, m for miu
   minimize(C'*x)
   subject to
     A'*L + m == C;
     A*x == b;
     diag(D) == x;
     D*m == t;
     x >= 0;
     m >= 0;
     x(1)+ x(2)+ x(3) -1 == 0;
cvx_end
end

%% Exercise 9.1.b
close all
clc
clear all

CT = [-1 -2 3];
C=transpose(CT);
A = [1 1 1];
b = 1;
D = [1 0 0; 0 1 0; 0 0 1]
t = 0;

cvx_begin sdp
   variables x(3) L(1,1) m(3,1)
   minimize(C'*x)
   subject to
     A'*L + m == C;
     A*x == b;
     diag(D) == x;
     D*m == 0;
     x >= 0;
     m >= 0;
     x(1)+ x(2)+ x(3) -1 == 0;
cvx_end

%% Exercise 9.1.c
close all
clc
clear all

CT = [-1 -1];
C=transpose(CT);
A = [1 -1];
b = 0;
D = [1 0; 0 1]
t = 0;

cvx_begin sdp
   variables x(2) L(1,1) m(2,1) % l for lambda, m for miu
   minimize(C'*x)
   subject to
     A'*L + m == C;
     A*x == b;
     diag(D) == x;
     D*m == 0;
     x >= 0;
     m >= 0;
     x(1)- x(2) == 0;
cvx_end
