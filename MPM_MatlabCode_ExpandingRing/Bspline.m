
function [phi, gphi, ctrlPs] = Bspline(x, deg, xp)

Nknots = length(x);
phi_vec = zeros(1,Nknots-1);
gphi_vec = zeros(1,Nknots-1);

kSpan = 0;
for i = 1:Nknots-1

    if (xp>=x(1,i) && xp<x(1,i+1))

        kSpan = i;
        phi_vec(1,kSpan) = 1;
        break;

    end

end

if kSpan == 0

    disp('**** ERRORRRR Unknown kSpan in Bspline() \n');

end

ctrlPs = kSpan-deg:kSpan;
set1 = [];

for p = 1:deg

    set1 = kSpan-deg:kSpan+deg-p;
    set2 = set1+1;
    set3 = set1+p;
    set4 = set2+p;
    Nset = size(set1,2);
    
    gphi_left = (p*ones(1,Nset))./(x(:,set3)-x(:,set1)).*phi_vec(:,set1);
    gphi_left(isnan(gphi_left)) = 0;
    phi_left = gphi_left.*(xp-x(:,set1))/p;
    gphi_right = (p*ones(1,Nset))./(x(:,set4)-x(:,set2)).*phi_vec(:,set2);
    gphi_right(isnan(gphi_right)) = 0;
    phi_right = gphi_right.*(x(:,set4)-xp)/p;
    
    phi_vec(1,set1) = phi_left+phi_right;
    gphi_vec(1,set1) = gphi_left-gphi_right;

end

phi = phi_vec(1,set1);
gphi = gphi_vec(1,set1);

end
