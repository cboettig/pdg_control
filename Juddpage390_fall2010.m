%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using the collocation method to solve boundary value problems
% This code solves the problem on page 390 in Judd, 
% which is a life-cycle consumption and savings model. 
%
% ARE 254 Dynamic Optimization Fall 2009 Updated Fall 2010
% Prof. J. Sanchirico
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Juddpage390_fall2010
clear all
close all

% Economic parameters
R=.1;
gam=-2;
rho1=.05;
W=1; % fixed wage

% Boundary conditions
A0=0; AF=0;

% Time span
t0=0; tf=50;

figure 
plot(1:tf,.5+([1:tf]/10)-4*([1:tf]/50).^2)
title('Wage profile over time')

% Chebyshev nodes
n=11;  % n-1 Degree of polynomial
m=n; % number of nodes 
a=t0; b=tf; 

for i=1:m
    xn(i)=cos(((m-i+.5)/m)*pi);
    Time(i)= (a+b)/2 + xn(i)*(b-a)/2;
end

% Generate the T(x) functions at the nodes with the coefficents
T=Chebybasis(n,m,xn);
dT=(2/(b-a))*Dchebybasis(n,m,xn,T);

c=zeros(n,2); c(1,:)=1; % initial guess
options=optimset('MaxFunEvals', 10e6, 'MaxIter', 1e4, 'LargeScale', 'on');%,...
[c,r1]=fsolve(@(c)RESID1(c,R,gam,rho1,W,A0,AF,n,a,b,xn,Time,T,dT),c',options);

Cc=T*c(1,:)'; 
Aa=T*c(2,:)'; 

% Extracting the elements of r1 that apply to state equations and not
% intitial conditions
R1_S=r1(2:m+1);
R2_S=r1(2+m:length(r1)-1);
Tplot=linspace(t0,tf,501);
Co=interp1(Time,Cc,Tplot);
As=interp1(Time,Aa,Tplot);

% Solving the residual equation at more points in the range [-1,1] using
% the collocation coefficents 
X=linspace(-1,1,500);
TIME=(a+b)/2 + X*(b-a)/2;
T=Chebybasis(n,length(X),X);
dT=(2/(b-a))*Dchebybasis(n,length(X),X,T);
[RES1 R1 R2]=RESID1(c,R,gam,rho1,W,A0,AF,n,a,b,X,TIME,T,dT);
    
figure
subplot(211)
plot(TIME, R1)
title('Residual on dc/dt: collocation')
subplot(212)
plot(TIME, R2)
title('Residual on dA/dt: collocation')


% plotting the solution
figure
subplot(211)
plot(Tplot,Co,'r')
title('Consumption with W(t): Collocation')
subplot(212)
plot(Tplot,As,'k')
title('Assets Solution with W(t): Collocation')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solving the W(t)= W using the built-in matlab solver that uses
% piecewise collocation methods (you can also solve W(t) varying by 
% uncommenting the equations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% using the chebyshev nodes for the mesh points
for i=1:n
    XX(i)=(a+b)/2 + ((b-a)/2)*cos(((n-i+.5)/n)*pi); % x is the time dimension
end

solinit = bvpinit(XX,[1 1]);
options = bvpset('Stats','on','RelTol',1e-6);
sol = bvp4c(@D1SD,@res1,solinit,options,R,gam,rho1,W,A0,AF);

x = sol.x; t=x;
y = sol.y; x=y';
figure
subplot(211)
plot(t,x(:,1),'r')
title('Consumption Solution with W(t)=W')
subplot(212)
plot(t,x(:,2),'k')
title('Asset Solution with W(t)=W')



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phase plane 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
C=linspace(0,5,20); 
A=linspace(-5,10,20);
plot(A,R*A+1,zeros(1,11),0:10,'g-','LineWidth',2)
hold on
% plotting the optimal solution for W(t) using collocation
plot(As,Co,'r--')
xlabel('A')
ylabel('C')
title(['C(0)=0 & C(T)=0 & T=' num2str(tf)])
% plotting the solution from BVP4c which uses W(t)=W.
plot(x(:,2),x(:,1))
legend('dA/dt', 'dc/dt','W(t)','W(t)=W')
% doing the directional arrows
[Ax,Cx]=meshgrid(A,C);
Axx=R*Ax-Cx+1;
Cxx=(1/gam)*Cx*(rho1-R);
quiver(Ax,Cx,Axx,Cxx)
axis([-5 10 0 5])

end


function [Res1 R1 R2]=RESID1(c,R,gam,rho1,W,A0,AF,n,a,b,xn,Time,T,dT)
 
co=zeros(1,length(c));
A=zeros(1,length(c));
 
co =  T*c(1,:)'; 
A =  T*c(2,:)';
cdot= dT*c(1,:)'; 
Adot= dT*c(2,:)';
 
R1= cdot - (1/gam)*co*(rho1-R);
R2= Adot - (R*A+(.5+(Time'./10)-4*(Time'/50).^2)-co); % Allowing Wage to
%vary
%R2= Adot - R*A - W + co; %fixed wage

Res1 = [A0-c(2,:)*((-1).^[1:length(c)])' R1' R2' AF-sum(c(2,:))];
 
end

function res = res1(ya,yb,R,gam,rho1,W,A0,AF)
res = [ya(2) - A0
       yb(2) - AF];
end

function dx = D1SD(t,x,R,gam,rho1,W,A0,AF)
% commented code lets you calculate the solution with W varying
% dx= [  (1/gam)*x(1)*(rho1-R)
%         R*x(2)-x(1)+(.5+(t./10)-4*(t./50).^2)];
   dx = [  (1/gam)*x(1)*(rho1-R)
            R*x(2)-x(1)+W];
end

function T=Chebybasis(n,m,x)
% Generate the T(x) Chebychev basis functions at the nodes
 for i=1:m
     for j=1:n
         T(i,j)=cos((j-1)*acos(x(i))); 
     end
 end
end

function dT=Dchebybasis(n,m,x,T)
% Generate the derivatives of the Cheby basis
dT=zeros(m,n);
 for i=1:length(x)
     for j=1:n
         dT(i,j)=(j-1)*sin((j-1)*acos(x(i)))/((1-x(i)^2)^.5);
     end
 end
end




