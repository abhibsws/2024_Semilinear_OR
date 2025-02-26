function uerr = DIRK_SL_PR_Eqn(TC,tf,s,p,q,scheme_no,dt,nt,lambda)
    % Parameter of the Prothero Robinson problem
    %lambda = -10^5;
    % DIRK scheme 
    [A,b,c] = SL_DIRK_Butcher(s,p,q,scheme_no);
    s = length(c); 
    tn = 0; y = utrue(TC,tn); % inital condition
    
    % time loop to calculate solutions
    for i=1:nt   
        t = tn+dt*c;
        g = zeros(1,s);     % for storing intermediate stages as column vectors
  	    R = zeros(1,s);                      % Evauation of rhs function at g_j
        for j = 1:s
            rhs = y-lambda*dt*A(j,j)*utrue(TC,t(j));
            if j>1
                rhs = rhs + dt*R(:,1:j-1)*(A(j,1:j-1))';
            end
            % solve for current stage by Newton iteration
            g_j = Newton_it(lambda,dt,A(j,j),rhs,y);
            % Calculating the right hand side
            R(:,j) = SL_PR_RHS(TC,t(j),g_j,lambda);
            g(:,j) = g_j;
        end
        if max(abs(A(end,:)-b))<1e-14
            y = g(:,s);
        else
            y = y + dt*(R*b');
        end
        tn = tn+dt;
    end
    uerr = abs(y-utrue(TC,tf));
end
      
%===========================Function routine==============================%
% Exact smooth solution
function z = utrue(TC,t)
    switch TC
        case 1, z =sqrt(1+t^2)-t; 
    end
end
% Right hand side of the problem
function z = SL_PR_RHS(TC,t,y,lambda)
    z = lambda*(y-utrue(TC,t))-2*y^2/(1+y^2);
end
% This is the nonlinear function that needs to be solved for stages
function z = nonlin_system_RHS(lambda,dt,ajj,g,rhs)
    z = (1-lambda*dt*ajj)*g + 2*dt*ajj*(g^2/(1+g^2)) - rhs;
end
% Newton method to solve stage values
function g = Newton_it(lambda,dt,ajj,rhs,u0)
    err = 1.0; tol = 1e-14;
    it = 1; maxit = 20;
    g = u0;
    while err>tol && it<maxit
        Pg = nonlin_system_RHS(lambda,dt,ajj,g,rhs);
        J = 1-lambda*dt*ajj + 4*dt*ajj*(g/(1+g^2)^2);
        dg = J\Pg;
        g = g - dg;
        err = max(abs(dg));
        %fprintf('Newton step %d, err = %10.6e \n',it,err);
        it = it+1;
    end
end
%===============================End=======================================%


