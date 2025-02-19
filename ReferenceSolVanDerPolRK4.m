function RefSol = ReferenceSolVanDerPolRK4(ep,dt,tf)
%-------------------------------------------------------------------------%
% This matlab function computes the numerical approximation for the problem
%------
%|x'(t) = y,
%|y'(t) = ep_inv*( (1-x^2)*y - x ).
%(x(0),y(0))=(2,-0.6666654321121172), ep_inv = 1e-5, and tf = 0.5,
%------
% by explicit RK4 method with dt = 10^(-6).
%-------------------------------------------------------------------------%
% Explicit RK4 Butcher tableau
A = [  0    0   0   0;
     1/2    0   0   0;
       0  1/2   0   0;
       0    0   1   0];
b = [1/6 1/3 1/3 1/6];
c = [0;1/2;1/2;1];

nt = ceil(tf/dt);
dt = tf/nt;
tn = 0;

ep_inv = ep^(-1);

% RefSol = zeros(2,nt+1);
RefSol(:,1) = [2;initial_cond(ep)]; % initial condition
u = RefSol(:,1);
% time loop to calculate solutions
for i=1:nt 
    t = tn+dt*c;
    g1 = VanDerPolrhs(t(1),u,ep_inv);
    g2 = VanDerPolrhs(t(2),(u+dt*0.5*g1),ep_inv);
    g3 = VanDerPolrhs(t(3),(u+dt*0.5*g2),ep_inv);
    g4 = VanDerPolrhs(t(4),(u+dt*g3),ep_inv);
    unew = u+dt*( b(1)*g1+b(2)*g2+b(3)*g3+b(4)*g4 );
    RefSol(:,i+1) = unew;
    tn = tn+dt;
    u = unew;
end
end

%===========================Function routine==============================%

function R = VanDerPolrhs(t,u,ep_inv)
    R(1,1) = u(2);
    R(2,1) = ep_inv*( (1-u(1)^2)*u(2) - u(1) );
end

function y0 = initial_cond(ep)
    y0 = -2/3 + (10/81)*ep - (292/2187)*ep^2 - (1814/19683)*ep^3;
end



