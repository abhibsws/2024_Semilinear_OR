xL = 0; xR = 1; m = 6;
[D,xold] = ChebArbDomain(m,xR,xL)

[w, x] = clenshawcurtiswx(m+1,xL,xR)

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

function [w,x] = clenshawcurtiswx(points, a, b)
% Computes Chebyshev (Clenshaw-Curtis) nodes and quadrature weights
% Input:  points - number of quadrature points (must be odd and ≥ 3)
% Output: x - Chebyshev nodes in [-1,1]
%         w - corresponding quadrature weights

    N = points - 1;
    
    n = (0:N/2)';
    k = 0:N/2;
    
    D = 2 * cos(2 * pi * (n * k) / N) / N;
    D(1, :) = 0.5 * D(1, :);  % halve the first row

    d = [1; (2 ./ (1 - (2:2:N).^2))'];
    
    w_half = D * d;
    w = [w_half; flipud(w_half(1:end-1))];
    
    x = cos((0:N)' * pi / N);

    % Scale points and weights to [a, b]
    x = (b - a) / 2 * x + (a + b) / 2;
    w = w * (b - a) / 2;
end