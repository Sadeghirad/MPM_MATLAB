
function b = Material(F,flag)

Lambda = 0;
mu = 5e6;

I=eye(2,2);
J=det(F);
LnJ=log(J);

if strcmp(flag,'Potential')
    
    b=0.5*Lambda*LnJ^2-mu*LnJ+0.5*mu*(trace(F'*F)-2);
    
elseif strcmp(flag,'Stress')
    
    b=Lambda/J*LnJ*I+mu/J*(F*F'-I);
    
end
