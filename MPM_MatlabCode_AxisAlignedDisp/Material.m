
function b = Material(F,flag)

Lambda = 5.769230769230768e+06;
mu = 3.846153846153846e+06;

I=eye(2,2);
J=det(F);
LnJ=log(J);

if strcmp(flag,'Potential')

    b=0.5*Lambda*LnJ^2-mu*LnJ+0.5*mu*(trace(F'*F)-2);

elseif strcmp(flag,'Stress')

    b=Lambda/J*LnJ*I+mu/J*(F*F'-I);

end
