function [uerr,duerr] = DIRKspqVisBurgerEqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt)
%=========================================================================%
% eqn = 'VisBurger'. 
% m = number of grid points without boundaries.
% spatial_order = is the central finite difference order to approximate
% spatial derivative.
% tf = final time.
% s = number of stages of DIRK scheme.
% p = order of DIRK scheme.
% q = weak stage order of DIRK scheme.
%=========================================================================%
    % Spatial approximation by finite difference 
    x = linspace(0,1,m+1)';
    xi = x(2:end-1);
    h = x(2)-x(1);
    e = ones(m-1,1);
    nu = 0.1;
    switch spatial_order
        case 2
            % M is for \nu*u_xx
            M = spdiags([e -2*e e],-1:1,m-1,m-1);
            M = nu*M ./h^2; M =sparse(M);
            % N is for u_x
            N = spdiags([-e 0*e e],-1:1,m-1,m-1);
            N = N./(2*h); N =sparse(N);
        case 4 
            % M is for \nu*u_xx
            M = spdiags([-e 16*e -30*e 16*e -e],-2:2,m-1,m-1);
            M(1,1:5)=[-15,-4,14,-6,1]; M(m-1,m-5:m-1) = [1,-6,14,-4,-15];
            M  = nu*M ./(12*h^2); M =sparse(M);
            % N is for u_x
            N = spdiags([e -8*e 0*e 8*e -e],-2:2,m-1,m-1);
            N(1,1:4) = [-10,18,-6,1]; N(m-1,m-4:m-1) = [-1,6,-18,10]; 
            N = N./(12*h); N =sparse(N);
        case 6
            % M is for \nu*u_xx
            M = spdiags([2*e -27*e 270*e -490*e 270*e -27*e 2*e],-3:3,m-1,m-1);
            M(1,1:7)= [-70 -486 855 -670 324 -90 11]; 
            M(m-1,m-7:m-1) = [11 -90 324 -670 855 -486 -70];
            M(2,1:7) = [214 -378 130 85 -54 16 -2];
            M(m-2,m-7:m-1) = [-2 16 -54 85 130 -378 214];
            M  = nu*M ./(180*h^2); M =sparse(M);
            % N is for u_x
            N = spdiags([-1*e 9*e -45*e 0*e 45*e -9*e 1*e],-3:3,m-1,m-1);
            N(1,1:6) = [-77,150,-100,50,-15,2];
            N(2,1:6) = [-24,-35,80,-30,8,-1];
            N(m-1,m-6:m-1) = [-2,15,-50,100,-150,77]; 
            N(m-2,m-6:m-1) = [1,-8,30,-80,35,24]; 
            N = N./(60*h); N =sparse(N);
    end
    
    % DIRK scheme 
    [A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no);
    %-------------------------------------------------------------------------%
    s = length(c);
    tn = 0.0;
    u = ue(TC,xi,tn); % Inital condition evaluation at tn = 0
    %Time loop till final time tf = 1
    for i = 1:nt
        t = tn+dt*c;
        g = zeros(m-1,s);     % for storing intermediate stages as column vectors
  	    R = zeros(m-1,s);                      % Evauation of rhs function at g_j
        for j = 1:s
            [bc_uxx,bc_ux] = G(TC,xi,t(j),h,spatial_order,nu);
            rhs = u + A(j,j)*dt*bc_uxx + A(j,j)*dt*F(TC,xi,t(j),nu); 
            if j>1
                rhs = rhs + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end
            % solve for the current stage by Newton iteration
            g_j = Newton_it(dt,A(j,j),M,N,bc_ux,rhs,u); 
            % Calculating MU(t)- U(t).*[NU(t) + bc_ux] + bc_uxx + F(t) at (tn+dt*c_j,g_j)
            R(:,j) = M*g_j - g_j.*( N*g_j + bc_ux ) + bc_uxx + F(TC,xi,t(j),nu);
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
    u_num = [ue(TC,x(1),tf);u;ue(TC,x(end),tf)]; du_num = duapprox(h,u_num,5);
    u_true = ue(TC,x,tf); du_true = ue_x(TC,x,tf);

    % max norm
    % uerr = max(abs(u_num-u_true)); duerr = max(abs(du_num-du_true));

    % l_2 norm
    uerr = sqrt(sum((u_num - u_true).^2) / (m+2));
    duerr = sqrt(sum((du_num - du_true).^2) / (m+2));
end

%------------------------Function Routine---------------------------------%
% Exact Solution
function z = ue(TC,x,t)
    switch TC
        case 1, C1 = 10; C2 = 20; z = cos(2+C1*t)*sin(0.2+C2*x);
        case 2, C1 = 5; C2 = 10; z = cos(2+C1*t)*sin(0.2+C2*x);
        case 3, C1 = 200; C2 = 50; z = cos(2+C1*t)*sin(0.2+C2*x);
        case 4, C1 = 1; z = cos(C1*t)*ones(size(x));
    end
end
% time derivative of exact solution
function z = ue_t(TC,x,t)
    switch TC
        case 1, C1 = 10; C2 = 20; z =  -C1*sin(2+C1*t)*sin(0.2+C2*x);
        case 2, C1 = 5; C2 = 10; z =  -C1*sin(2+C1*t)*sin(0.2+C2*x);
        case 3, C1 = 200; C2 = 50; z =  -C1*sin(2+C1*t)*sin(0.2+C2*x);
        case 4, C1 = 1; z = -C1*sin(C1*t)*ones(size(x));
    end
end
%First order spatial derivative of exact solution
function z = ue_x(TC,x,t)
    switch TC
        case 1, C1 = 10; C2 = 20; z =  cos(2+C1*t)*C2*cos(0.2+C2*x);
        case 2, C1 = 5; C2 = 10; z =  cos(2+C1*t)*C2*cos(0.2+C2*x);
        case 3, C1 = 200; C2 = 50; z =  cos(2+C1*t)*C2*cos(0.2+C2*x);
        case 4, C1 = 1; z = cos(C1*t)*zeros(size(x)); 
    end
end
%Second order spatial derivative of exact solution
function z = ue_xx(TC,x,t)
    switch TC
        case 1, C1 = 10; C2 = 20; z =  -cos(2+C1*t)*C2^2*sin(0.2+C2*x);
        case 2, C1 = 5; C2 = 10; z =  -cos(2+C1*t)*C2^2*sin(0.2+C2*x);
        case 3, C1 = 200; C2 = 50; z =  -cos(2+C1*t)*C2^2*sin(0.2+C2*x);
        case 4, C1 = 1; z = cos(C1*t)*zeros(size(x));
    end
end
% Forcing function f(x,t) = u_t+uu_x-\nu*u_xx as vector
function z = F(TC,x,t,nu)
    z = ue_t(TC,x,t)+ue(TC,x,t).*ue_x(TC,x,t)-nu*ue_xx(TC,x,t);
end
% Exact boundary ondition as vectors
function [bc_uxx,bc_ux] = G(TC,xi,t,h,spatial_order,nu)
    bc_uxx = zeros(length(xi),1);
    bc_ux = zeros(length(xi),1);
    switch spatial_order
        case 2
            bc_uxx(1,1) = nu*g0(TC,t)/h^2; bc_uxx(end,1) = nu*g1(TC,t)/h^2;
            
            bc_ux(1,1) = -1*g0(TC,t)/(2*h); bc_ux(end,1) = 1*g1(TC,t)/(2*h);
        case 4
            bc_uxx(1,1)= nu*10*g0(TC,t)/(12*h^2); bc_uxx(end,1) = nu*10*g1(TC,t)/(12*h^2);
            bc_uxx(2,1)= nu*-1*g0(TC,t)/(12*h^2); bc_uxx(end-1,1) = nu*-1*g1(TC,t)/(12*h^2);
             
            bc_ux(1,1) = g0(TC,t)*(-3/(12*h)); bc_ux(end,1)  = g1(TC,t)*(3/(12*h));
            bc_ux(2,1) = g0(TC,t)*(1/(12*h));  bc_ux(end-1,1)  = g1(TC,t)*(-1/(12*h));
        case 6
            bc_uxx(1,1)= nu*126*g0(TC,t)/(180*h^2); bc_uxx(end,1) = nu*126*g1(TC,t)/(180*h^2);
            bc_uxx(2,1)= nu*-11*g0(TC,t)/(180*h^2); bc_uxx(end-1,1) = nu*-11*g1(TC,t)/(180*h^2);
            bc_uxx(3,1)= nu*2*g0(TC,t)/(180*h^2); bc_uxx(end-2,1) = nu*2*g1(TC,t)/(180*h^2);
            
            bc_ux(1,1)= -10*g0(TC,t)/(60*h); bc_ux(end,1) = 10*g1(TC,t)/(60*h);
            bc_ux(2,1)= 2*g0(TC,t)/(60*h); bc_ux(end-1,1) = -2*g1(TC,t)/(60*h);
            bc_ux(3,1)= -1*g0(TC,t)/(60*h); bc_ux(end-2,1) = 1*g1(TC,t)/(60*h);
    end
end
% Exact boundary condition left
function z = g0(TC,t)
    z = ue(TC,0,t);
end
% Exact boundary condition right
function z = g1(TC,t)
    z = ue(TC,1,t);
end
% Newton iteration to solve stage involving nonlinearity
function w = Newton_it(dt,ajj,M,N,bc_ux,rhs,u0)
    n = size(M,2); I = speye(n,n);
    err = 1.0; tol = 1e-12;
    it = 1; maxit = 10;
    w = u0;
    while err>tol && it<maxit
        Pw = (I-dt*ajj*M)*w + dt*ajj*w.*(N*w + bc_ux ) - rhs;
        J = (I-dt*ajj*M) + dt*ajj*( w.*N + diag(N*w+bc_ux) ); J = sparse(J);
        dw = J\Pw;
        w = w-dw;
        err = max(abs(dw));
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
%=========================================================================%