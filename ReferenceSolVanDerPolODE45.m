function RefSol = ReferenceSolVanDerPolODE45(ep,tf)
%-------------------------------------------------------------------------%
% This matlab code solves the Van der Pol oscillator with ode45
%------
%|x'(t) = y,
%|y'(t) = ep_inv*( (1-x^2)*y - x ).
%(x(0),y(0))=(2,-2/3 + (10/81)*epsilon - (292/2187)*epsilon^2 - (1814/19683)*epsilon^3) 
% and tf = 0.5.
%------
%-------------------------------------------------------------------------%
    % Time span
    tspan = [0, tf];
    
    % Initial conditions
    y1_0 = 2;
    y2_0 = -2/3 + (10/81)*ep - (292/2187)*ep^2 - (1814/19683)*ep^3;
    y0 = [y1_0; y2_0];
    
    % ODE system
    vdp = @(t, y) [...
        y(2);
        (1/ep)*((1 - y(1)^2)*y(2) - y(1))
    ];
    
    % Solver options: tight tolerance
    options = odeset('RelTol',5e-14,'AbsTol',5e-14);
    
    % Solve
    [t, y] = ode45(vdp, tspan, y0, options);
    
    RefSol = y(end,:)';
end