function uerr = DIRKspqVanDerPol(TC,tf,s,p,q,scheme_no,dt,nt,ep)
%-------------------------------------------------------------------------%
% This matlab function computes the numerical approximation for the problem
%------
%|x'(t) = y,
%|y'(t) = ep^{-1}*( (1-x^2)*y - x ).
%(x(0),y(0))=(2,-0.6666654321121172), and tf = 0.5,
%------
% using by a choice of s-stages pth order scheme with WSO q. 
% ---Input---
    % TC = choice of initial condition case
    % tf = final time
    % (s,p,q) = (stage, order, wso) for DIRK scheme
    % scheme numner = scheme_no
    % dt = time step size and nt = number of time steps
% ---Output---
    % Error at the final time as uerr.
%-------------------------------------------------------------------------%
% Parameter for Van Der Pol equation
%ep = 1e-5; 
ep_inv = ep^(-1);
% DIRK scheme 
[A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no);
tn = 0;
switch TC
    case 1       
        u = [2;initial_cond(ep)]; % initial condition
end
sol = zeros(2,nt+1);
sol(:,1) = u;
% time loop to calculate solutions
for i=1:nt   
    t = tn+dt*c;
    g = zeros(2,s);     % for storing intermediate stages as column vectors
  	R = zeros(2,s);                      % Evauation of rhs function at g_j
    for j = 1:s
        rhs = u;
        if j>1
            rhs = rhs + dt*R(:,1:j-1)*(A(j,1:j-1))';
        end
        % solve for current stage by Newton iteration
        g_j = Newton_it(t(j),dt,ep_inv,A(j,j),rhs,u);
        % Calculating the right hand side
        R(:,j) = VanDerPolrhs(t(j),g_j,ep_inv);
        g(:,j) = g_j;
    end
    if max(abs(A(end,:)-b))<1e-14
        u = g(:,s);
    else
        u = u + dt*(R*b');
    end
    tn = tn+dt;
    sol(:,i+1) = u;
end
% true_sol is copied from RefSol...(:,end)
if ep == 1e-1
    %true_sol = [1.613276839978073; -0.943670141852980]; % ep = 1e-1 and dt = 1e-7
     true_sol = [1.613276839978089; -0.943670141852946]; % by ODE45 with tol 5e-14
elseif ep == 1e-2
    %true_sol = [1.598829069390465; -1.018139709121219]; % ep = 1e-2 and dt = 1e-7
     true_sol = [1.598829069390726; -1.018139709120860]; % by ODE45 with tol 5e-14
elseif ep == 1e-3
    %true_sol = [1.596980778659228; -1.029103015879397]; % ep = 1e-3 and dt = 1e-7
     true_sol = [1.596980778659659; -1.029103015878781]; % by ODE45 with tol 5e-14
elseif ep == 1e-4
    %true_sol = [1.596789700158140; -1.030263287387102]; % ep = 1e-4 and dt = 1e-7
     true_sol = [1.596789700158133; -1.030263287387118]; % by ODE45 with tol 5e-14
elseif ep == 1e-5
    %true_sol = [1.596770525704857; -1.030380015613961]; % ep = 1e-5 and dt = 1e-7
     true_sol = [1.596770525704804; -1.030380015614045 ]; % by ODE45 with tol 5e-14
elseif ep == 1e-6
    %true_sol = [1.596768607589027; -1.030391695517091]; % ep = 1e-6 and dt = 1e-7
     true_sol = [1.596768607588914; -1.030391695517268]; % by ODE45 with tol 5e-14
     
else
    fprintf('The choice of the epsilon is NOT coded. \n')
end

% max norm
%uerr = max(abs(true_sol-sol(:,end)));

%l2 norm
uerr = sqrt(sum((true_sol - sol(:, end)).^2) / numel(true_sol));
%----------------------------------------------
% This step loads file every time so we do it another way above
% solfile = sprintf('RefSol_mu%d_tf%d',mu,tf);
% load(solfile)
% uerr = max(abs(RefSol_mu500_tf10(:,end)-sol(:,end)));
%---------------------------------------------
end    

%===========================Function routine==============================%

function R = VanDerPolrhs(t,u,ep_inv)
    R(1,1) = u(2,1);
    R(2,1) = ep_inv*( (1-u(1,1)^2)*u(2,1) - u(1,1) );
end

function P = nonlin_system(t,g,dt,ep_inv,ajj,rhs)
    P(1,1) = g(1,1)-dt*ajj*g(2,1)-rhs(1,1);
    P(2,1) = g(2,1)-ep_inv*dt*ajj*(1-g(1,1)^2)*g(2,1)+ep_inv*dt*ajj*g(1,1)-rhs(2,1);
end

function y0 = initial_cond(ep)
    y0 = -2/3 + (10/81)*ep - (292/2187)*ep^2 - (1814/19683)*ep^3;
end


function g = Newton_it(t,dt,ep_inv,ajj,rhs,u0)
    err = 1.0; tol = 1e-12;
    it = 1; maxit = 10;
    g = u0;
    while err>tol && it<maxit
        Pg = nonlin_system(t,g,dt,ep_inv,ajj,rhs);
        J = [               1                            -dt*ajj;
              ep_inv*dt*ajj*(1+2*g(1,1)*g(2,1))     1-ep_inv*dt*ajj*(1-g(1,1)^2)];
        dg = J\Pg;
        g = g - dg;
        err = max(abs(dg));
%         fprintf('Newton step %d, err = %10.6e \n',it,err);
        it = it+1;
    end
end


