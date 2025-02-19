function [uerr,duerr] = DIRK_SemiLinAdvEqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt)
    % Spatial approximation by finite difference 
    x = linspace(0,1,m+1)';
    xi = x(2:m+1);
    h = x(2)-x(1);
    e = ones(m,1);
    switch spatial_order
        case 2
            % N is for u_x
            N = spdiags([-e 0*e e],-1:1,m,m);
            N(m,m-2:m)=[1 -4 3];
            N = N./(2*h); N =sparse(N);
        case 4 
            % N is for u_x
            N = spdiags([e -8*e 0*e 8*e -e],-2:2,m,m);     
            N(1,1:4) = [-10,18,-6,1];
            N(m-1,m-4:m) = [-1,6,-18,10,3];
            N(m,m-4:m) = [3,-16,36,-48,25];
            N = N./(12*h); N =sparse(N);  
        case 5
             %--------------------------------%
             % N is for u_x: 6th order in the interior and 6th order near bdd
             % will produce unstable spatial discretization. Similarly 6th
             % order in the interior and 5th order near bdd will also lead to
             % unstable discretization. But 6th order in the interior and 4th 
             % order near bdd produces an overall 5th orderstable spatial 
             % discretization and we will use this semi-discretization for our
             % high order time-stepping methods. For the details experiments
             % look at the overleaf document and the code from 
             % Current Projects/ExplicitWSO/NumericalResultsForPaper.
            %---------------------------------%
            %N is for u_x: 6th order in the interior and 4th order near bdd
            N = spdiags([-1*e 9*e -45*e 0*e 45*e -9*e 1*e],-3:3,m,m);
            N(1,1:5) = [-50,    90,   -30,     5, 0]; % x = x1, 4th order one-sided FD is used (require 5 point stencil)
            N(2,1:5) = [-40,     0,    40,    -5, 0];  % x = x2, 4th order one-sided FD is used
            N(m-2,m-5:m) = [0, 5,   -40,     0,    40,    -5]; % x = x_{m-2}, 4th order one-sided FD is used
            N(m-1,m-4:m) = [-5,    30,   -90,    50,    15]; % x = x_{m-1}, 4th order one-sided FD is used
            N(m,m-4:m) = [15,   -80,   180,  -240,   125];   % x = x_{m} = 1, 4th order one-sided FD is used    
            N = N./(60*h); N =sparse(N);
    end
    
    % DIRK scheme 
    [A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no);
    %--------------------------------------------------------------------------
    s = length(c); tn = 0.0; 
    u = ue(TC,xi,tn); % Inital condition evaluation at tn = 0
    % Time loop
    for i = 1:nt
        t = tn+dt*c;
        g = zeros(m,s);     % for storing intermediate stages as column vectors
        R = zeros(m,s);                      % Evauation of rhs function at g_j
        for j = 1:s
            bc_ux = G(TC,xi,t(j),h,spatial_order);
            rhs = u - A(j,j)*dt*bc_ux; 
            if j>1
                rhs = rhs + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end
            % solve for current stage by Newton iteration
            g_j = Newton_it(dt,A(j,j),N,rhs,u);
            % Calculating  -( NU(t) + bc_ux ) + F(U(t))
            R(:,j) = -( N*g_j + bc_ux ) + g_j.^2;
            g(:,j) = g_j;
        end
        if max(abs(A(end,:)-b))<1e-14
            u = g(:,s);
        else
            u = u + dt*(R*b');
        end
        tn = tn+dt;
    end
    % Computing error
    u_num = [ue(TC,x(1),tf);u]; du_num = duapprox(h,u_num,5);
    u_true = ue(TC,x,tf); du_true = ue_x(TC,x,tf);

    % max norm
    %uerr = max(abs(u_num-u_true)); duerr = max(abs(du_num-du_true));

    % l_2 norm
    uerr = sqrt(sum((u_num - u_true).^2) / (m+1));
    duerr = sqrt(sum((du_num - du_true).^2) / (m+1));
end

%------------------------Function Routine---------------------------------%
% Exact Solution
function z = ue(TC,x,t)
    switch TC
        case 1, z = (sin(pi*(x-t)).*sin(pi*(x-t)))./...
                (1-t*sin(pi*(x-t)).*sin(pi*(x-t)));
    end
end
% Derivative of the Exact Solution
function z = ue_x(TC,x,t)
    switch TC
        case 1, z = -(2*pi*cos(pi*(t - x)).*sin(pi*(t - x)))./(t*sin(pi*(t - x)).^2 - 1).^2;
    end
end
% Exact boundary ondition as vectors
function bc_ux = G(TC,xi,t,h,spatial_order)
    bc_ux = zeros(length(xi),1);
    switch spatial_order
        case 2
          bc_ux(1,1) = -1*g0(TC,t)/(2*h);
        case 4
            bc_ux(1,1) = g0(TC,t)*(-3/(12*h)); bc_ux(2,1) = g0(TC,t)*(1/(12*h)); 
        case 5
           %6th order in the interior and 4th order near bdd
            bc_ux(1,1) = g0(TC,t)*(-15/(60*h)); bc_ux(2,1) = g0(TC,t)*(5/(60*h)); bc_ux(3,1) = g0(TC,t)*(-1/(60*h)); 
    end
end
% Exact boundary condition left (only need this Dirchlet bc as the information
% is moving to the right  with velocity 1).
function z = g0(TC,t)
    z = ue(TC,0,t);
end
% Newton iteration to solve stage involving nonlinearity
function g = Newton_it(dt,ajj,N,rhs,u0)
    m = size(N,2); I = speye(m,m);
    err = 1.0; tol = 1e-14;
    it = 1; maxit = 10;
    g = u0;
    while err>tol && it<maxit
        Pg = (I+dt*ajj*N)*g - dt*ajj*g.^2 - rhs;
        J = (I+dt*ajj*N) - 2*dt*ajj*diag(g); J = sparse(J);
        dg = J\Pg;
        g = g-dg;
        err = max(abs(dg));
        %fprintf('newton it. step %d, err = %10.6e\n',it,err);
        it = it+1;
    end
end
% Different order derivative approximation from function values
function du = duapprox(dx,u,order)
   % compute u_x from u
    switch order 
        case 2 % 2nd order approx.
            du = [sum([-3;4;-1].*u(1:3))./(2*dx);...
                  (u(3:end)-u(1:end-2))./(2*dx);...
                  sum([1;-4;3].*u(end-2:end))./(2*dx)];

        case 4 % 4th order approx.
            du = [sum([-25;48;-36;16;-3].*u(1:5))./(12*dx);...
                  sum([-3;-10;18;-6;1].*u(1:5))./(12*dx);...
                  (u(1:end-4)-8*u(2:end-3)+8*u(4:end-1)-u(5:end))./(12*dx);...
                  sum([-1;6;-18;10;3].*u(end-4:end))./(12*dx);...
                  sum([3;-16;36;-48;25].*u(end-4:end))./(12*dx)];
              
        case 5 % 6th order approx.
            du = [sum([-147;360;-450;400;-225;72;-10].*u(1:7))./(60*dx);...
                  sum([-10;-77;150;-100;50;-15;2].*u(1:7))./(60*dx);...
                  sum([2;-24;-35;80;-30;8;-1].*u(1:7))./(60*dx);...
                  (-1*u(1:end-6)+9*u(2:end-5)-45*u(3:end-4)+45*u(5:end-2)-9*u(6:end-1)+1*u(7:end))./(60*dx);...
                  sum([1;-8;30;-80;35;24;-2].*u(end-6:end))./(60*dx);...
                  sum([-2;15;-50;100;-150;77;10].*u(end-6:end))./(60*dx);...
                  sum([10;-72;225;-400;450;-360;147].*u(end-6:end))./(60*dx)];
    end
end
%-------------------------------------------------------------------------%