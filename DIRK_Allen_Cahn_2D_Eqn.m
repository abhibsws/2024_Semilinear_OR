function err = DIRK_Allen_Cahn_2D_Eqn(TC,tf,m,spatial_order,s,p,q,scheme_no,dt,nt)
    % parameters
    alp = 0.1; beta = 3; 
 
    % Space
    xL = 0; xR = 1; yL = xL; yR = xR;
    switch spatial_order
        case 'Cheb'
            [D,x] = ChebArbDomain(m,xR,xL); y = x;
    end
    
    % Laplacian
    D2 = D^2; I1D = eye(m+1);
    M = kron(I1D,D2) + kron(D2,I1D); M = alp*sparse(M);         
    
    [X,Y] = meshgrid(x,y); xx = X(:); yy = Y(:);
    
    B = find(abs(xx)==xL | abs(xx)==xR | abs(yy)==yL | abs(yy)==yR);% boundary points' indexes
    
    
    % DIRK scheme 
    [A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no); 
    
    % Time stepping
    s = length(c); 
    tn = 0.0; 
    u = ue(TC,xx,yy,tn);
    
    for i = 1:nt
        t = tn+dt*c;
        g = zeros((m+1)^2,s); % for storing intermediate stages as column vectors
        R = zeros((m+1)^2,s);     % Evauation of non-stiff non-linear part at g_j
        for j = 1:s
            f_val = F(TC,alp,beta,xx,yy,t(j));
            rhs = u+dt*A(j,j)*f_val; 
            if j>1
                rhs = rhs + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end 
            % Solve for the current stage using Newton iteration: the boundary condition
            % is incorporated within the nonlinear solver. The initial guess contains
            % the information of the Dirichlet boundary conditions, which are imposed
            % by replacing the indices corresponding to the boundary points with
            % the exact values of the solution at time t(j) for the current stage. These
            % values must be carried over throughout all calculations within the nonlinear
            % solver during the Newton iteration. Therefore, when we solve for the stages
            % (w = w - dw) using the Jacobian matrix, we need to solve for the increment dw.
            % To solve for the increment, we replace the rows indexed by B in the Jacobian
            % matrix with the corresponding rows from an identity matrix of the same size and
            % set the corresponding elements of the right-hand side to zero, ensuring that
            % we get a zero solution for the indices B in dw. Consequently, w = w - dw will
            % preserve the Dirichlet boundary conditions.
            u0 = u;
            u0(B,1) = ue(TC,xx(B),yy(B),t(j)); 
            g_j = Newton_it(B,dt,A(j,j),beta,M,rhs,u0); 
            % Calculating RHS of the system
            R(:,j) = M*g_j + nonlin_fun(beta,g_j)+F(TC,alp,beta,xx,yy,t(j));
            g(:,j) = g_j;
        end
        if max(abs(A(end,:)-b))<1e-14
            u = g(:,s);
        else
            u = u + dt*(R*b');
        end
        tn = tn+dt;
    end
    ut =  ue(TC,xx,yy,tf);
    % max norm
    % err = max(max(abs(u-ut)));

    % l2 norm
    err = sqrt(sum(sum((u - ut).^2)) / (m+1)^2);

    % Plot
    % ax = [0,1,0,1,[1,3]]; % plot axes
    % clf
    % subplot(1,2,1)
    % mesh(X,Y,reshape(ut,m+1,m+1)), axis(ax)
    % title(sprintf('True solution at t=%0.3f',tn))
    % subplot(1,2,2)
    % mesh(X,Y,reshape(u,m+1,m+1)), axis(ax)
end


%--------------------------------Functions--------------------------------%
% Exact Solution
function u = ue(TC,x,y,t)
    switch TC
        case 1, u = 2 + sin(2*pi*(x-t)).*cos(3*pi*(y-t));
    end
end

function u = ue_t(TC,x,y,t)
    switch TC
        case 1, u = 3*pi*sin(2*pi*(t - x)).*sin(3*pi*(t - y)) - 2*pi*cos(2*pi*(t - x)).*cos(3*pi*(t - y));
    end
end

function u = Lap_u(TC,x,y,t)
    switch TC
        case 1, u = 13*pi^2*cos(3*pi*(t - y)).*sin(2*pi*(t - x));
    end
end

% Nonlinear function
function z = nonlin_fun(beta,u)
    % u is expected to be a long vector of size (m+1)^2 by 1.
    z = beta*(u - u.^3);
end

% Forcing function f(x,y,t) = u_t-alp**(u_xx+u_yy)-beta*(u-u^3) 
function z = F(TC,alp,beta,x,y,t)
    % u is expected to be a long vector of size (m+1)^2 by 1.
    z = ue_t(TC,x,y,t)-alp*Lap_u(TC,x,y,t)-beta*(ue(TC,x,y,t)-ue(TC,x,y,t).^3);
end

% CHEB  compute D = differentiation matrix, x = Chebyshev grid
function [D,x] = cheb(N)
  if N==0, D=0; x=1; return, end
  x = cos(pi*(0:N)/N)'; 
  c = [2; ones(N-1,1); 2].*(-1).^(0:N)';
  X = repmat(x,1,N+1);
  dX = X-X';                  
  D  = (c*(1./c)')./(dX+(eye(N+1)));                 % off-diagonal entries
  D  = D - diag(sum(D'));                                % diagonal entries
end
  
function [D,x] = ChebArbDomain(N,xR,xL)
% when we solve on an interval [xL,xR] other than [-1,1]. Notice the order
% of the interval and how the end point is taken first and starting point
% at the end. 
  [Dcheb,scheb] = cheb(N);
  k1 = (xR - xL)/2;
  k2 = (xR + xL)/2;
  x = k2 + k1*scheb;
  D = Dcheb/k1;
end

% Newton iteration to solve stage involving nonlinearity
function w = Newton_it(B,dt,ajj,beta,M,rhs,u0)
    n = size(M,2); I = speye(n,n);
    err = 1.0; tol = 1e-12;
    it = 1; maxit = 10;
    w = u0;
    while err>tol && it<maxit
        Pw = (I-dt*ajj*M)*w - dt*ajj*nonlin_fun(beta,w) - rhs;
        Pw(B,1) = zeros(length(B),1); % modification for BCs
        J = (I-dt*ajj*M) - dt*ajj*beta*(I-3*diag(w.^2));
        J(B,:) = I(B,:); J = sparse(J);  % modification for BCs
        dw = J\Pw;
        w = w-dw;
        err = max(abs(dw));
        %fprintf('newton it. step %d, err = %10.6e\n',it,err);
        it = it+1;
    end
end
%-------------------------------------------------------------------------%