function b = AnalyticalSolution(X, Y, t, flag)

global interpolator gridSpacing

if strcmp(interpolator,'CPDI')
    
    X=X-gridSpacing; Y=Y-gridSpacing;
    
end

A = 0.05; % displacement amplitude
E = 1e7;  % modulus of elasticity
nu = 0.3; % Poisson's ratio
p = 1000; % density
C = 100;  % sqrt(E/rho)
Lambda = E/(1+nu)/(1-2*nu)*nu;
mu = E/2/(1+nu);

if strcmp(flag,'BodyForce')
    
    Fxx=1+2*A*pi*cos(2*pi*X)*sin(C*pi*t);
    Fyy=1+2*A*pi*cos(2*pi*Y)*sin(pi+C*pi*t);
    K=log(Fxx*Fyy)-1;
    u=[A*sin(2*pi*X)*sin(C*pi*t);A*sin(2*pi*Y)*sin(pi+C*pi*t)];
    b=pi^2*[u(1)*(4*mu/p-C^2-4*(Lambda*K-mu)/p/Fxx^2);u(2)*(4*mu/p-C^2-4*(Lambda*K-mu)/p/Fyy^2)];
    
elseif strcmp(flag,'Displacement')
    
    b=[A*sin(2*pi*X)*sin(C*pi*t);A*sin(2*pi*Y)*sin(pi+C*pi*t)];
    
elseif strcmp(flag,'DeformationGrad')
    
    Fxx=1+2*A*pi*cos(2*pi*X)*sin(C*pi*t);
    Fyy=1+2*A*pi*cos(2*pi*Y)*sin(pi+C*pi*t);
    b=[Fxx 0;0 Fyy];
    
elseif strcmp(flag,'Velocity')
    
    b=[A*sin(2*pi*X)*C*pi*cos(C*pi*t);A*sin(2*pi*Y)*C*pi*cos(pi+C*pi*t)];
    
end
