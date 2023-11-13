
function [nn,nodes,phi,nng,gnodes,gphi] = CalcSFGradSF(x,nCelly,conn,xp,F,InitParVec1,InitParVec2,deg,knotx,knoty,NctrlPx,snode)

global interpolator gridSpacing

numpar = size(xp,2);
nodes = zeros(snode,numpar);
phi = zeros(numpar,snode);
gnodes = zeros(snode,numpar);
gphi = zeros(2,snode,numpar);
nn = zeros(numpar,1);
nng = zeros(numpar,1);

if strcmpi(interpolator,'CPDI')

    for i = 1:numpar

        [nodes(:,i),phi(i,:)] = CalcSF_CPDI(xp(:,i),F(:,:,i),InitParVec1(:,i),InitParVec2(:,i),gridSpacing,x,conn,nCelly);
        nn(i) = numel(find(nodes(:,i)~=0));
        [gnodes(:,i),gphi(:,:,i)] = CalcGradSF_CPDI(xp(:,i),F(:,:,i),InitParVec1(:,i),InitParVec2(:,i),gridSpacing,x,conn,nCelly);
        nng(i) = numel(find(gnodes(:,i)~=0));

    end

elseif strcmpi(interpolator,'BSMPM')

    ctrlPx = zeros(numpar,deg+1);
    N_vecx = zeros(numpar,deg+1);
    B_vecx = zeros(numpar,deg+1);
    ctrlPy = ctrlPx;
    N_vecy = N_vecx;
    B_vecy = B_vecx;

    for i = 1:numpar

        [N_vecx(i,:),B_vecx(i,:),ctrlPx(i,:)] = Bspline(knotx,deg,xp(1,i)');
        [N_vecy(i,:),B_vecy(i,:),ctrlPy(i,:)] = Bspline(knoty,deg,xp(2,i)');

    end

    for i = 1:numpar

        nodes(:,i) = bspline_ctrlP_ids(ctrlPx(i,:),ctrlPy(i,:),NctrlPx,deg);
        phi(i,:) = bspline_ctrlP_sfs(N_vecx(i,:),N_vecy(i,:),deg);
        gphi(:,:,i) = bspline_ctrlP_dsfs(N_vecx(i,:),N_vecy(i,:),B_vecx(i,:),B_vecy(i,:),deg);
        nn(i) = snode;

    end

    nng = nn;
    gnodes = nodes;

elseif strcmpi(interpolator,'BSCPDI')

    for i = 1:numpar

        ParVec1 = F(:,:,i)*InitParVec1(:,i);
        ParVec2 = F(:,:,i)*InitParVec2(:,i);
        [xEx1,xEx2,xEx3,xEx4,ParVec1,ParVec2] = bsParDom(ParVec1,ParVec2,knotx,knoty,xp(:,i));
        [phi1,nodes1] = bsCorner(xEx1,deg,knotx,knoty,NctrlPx);
        [phi2,nodes2] = bsCorner(xEx2,deg,knotx,knoty,NctrlPx);
        [phi3,nodes3] = bsCorner(xEx3,deg,knotx,knoty,NctrlPx);
        [phi4,nodes4] = bsCorner(xEx4,deg,knotx,knoty,NctrlPx);
        vol1 = 4.*(ParVec1(1)*ParVec2(2)-ParVec1(2)*ParVec2(1));
        phi(i,:) = [phi1 phi2 phi3 phi4]/4;
        gphi(:,:,i) = [phi1*(ParVec1(2)-ParVec2(2))  phi2*(ParVec1(2)+ParVec2(2)) ...
                      -phi3*(ParVec1(2)-ParVec2(2)) -phi4*(ParVec1(2)+ParVec2(2))
                       phi1*(ParVec2(1)-ParVec1(1))  phi2*(-ParVec1(1)-ParVec2(1)) ...
                      -phi3*(ParVec2(1)-ParVec1(1)) -phi4*(-ParVec1(1)-ParVec2(1))]/vol1;
        nodes(:,i) = [nodes1;nodes2;nodes3;nodes4];
        nn(i) = snode;

    end

    nng = nn;
    gnodes = nodes;

end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  Auxiliary functions  %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gridNodes,phi] = CalcSF_CPDI(xp,F,InitParVec1,InitParVec2,gridSpacing,x,conn,nCelly)

    ParVec1 = F*InitParVec1;
    ParVec2 = F*InitParVec2;

    xEx = zeros(2,4);
    xEx(:,1) = [xp(1)-ParVec1(1)-ParVec2(1);xp(2)-ParVec1(2)-ParVec2(2)];
    xEx(:,2) = [xp(1)+ParVec1(1)-ParVec2(1);xp(2)+ParVec1(2)-ParVec2(2)];
    xEx(:,3) = [xp(1)+ParVec1(1)+ParVec2(1);xp(2)+ParVec1(2)+ParVec2(2)];
    xEx(:,4) = [xp(1)-ParVec1(1)+ParVec2(1);xp(2)-ParVec1(2)+ParVec2(2)];

    ie1 = (ceil(xEx(1,1)/gridSpacing))*nCelly-floor(xEx(2,1)/gridSpacing);
    ie2 = (ceil(xEx(1,2)/gridSpacing))*nCelly-floor(xEx(2,2)/gridSpacing);
    ie3 = (ceil(xEx(1,3)/gridSpacing))*nCelly-floor(xEx(2,3)/gridSpacing);
    ie4 = (ceil(xEx(1,4)/gridSpacing))*nCelly-floor(xEx(2,4)/gridSpacing);
    
    gridNodes = [conn(:,ie1);conn(:,ie2);conn(:,ie3);conn(:,ie4)];

    xpe1 = [2*(xEx(1,1)-x(1,conn(1,ie1)))/gridSpacing-1;2*(xEx(2,1)-x(2,conn(1,ie1)))/gridSpacing-1];
    xpe2 = [2*(xEx(1,2)-x(1,conn(1,ie2)))/gridSpacing-1;2*(xEx(2,2)-x(2,conn(1,ie2)))/gridSpacing-1];
    xpe3 = [2*(xEx(1,3)-x(1,conn(1,ie3)))/gridSpacing-1;2*(xEx(2,3)-x(2,conn(1,ie3)))/gridSpacing-1];
    xpe4 = [2*(xEx(1,4)-x(1,conn(1,ie4)))/gridSpacing-1;2*(xEx(2,4)-x(2,conn(1,ie4)))/gridSpacing-1];

    phis1 = 1/4*[(1-xpe1(1))*(1-xpe1(2)) (1+xpe1(1))*(1-xpe1(2)) ...
                 (1+xpe1(1))*(1+xpe1(2)) (1-xpe1(1))*(1+xpe1(2))];
    phis2 = 1/4*[(1-xpe2(1))*(1-xpe2(2)) (1+xpe2(1))*(1-xpe2(2)) ...
                 (1+xpe2(1))*(1+xpe2(2)) (1-xpe2(1))*(1+xpe2(2))];
    phis3 = 1/4*[(1-xpe3(1))*(1-xpe3(2)) (1+xpe3(1))*(1-xpe3(2)) ...
                 (1+xpe3(1))*(1+xpe3(2)) (1-xpe3(1))*(1+xpe3(2))];
    phis4 = 1/4*[(1-xpe4(1))*(1-xpe4(2)) (1+xpe4(1))*(1-xpe4(2)) ...
                 (1+xpe4(1))*(1+xpe4(2)) (1-xpe4(1))*(1+xpe4(2))];

    phi = [phis1 phis2 phis3 phis4]/4;
    
end


function [gridNodes,gphi] = CalcGradSF_CPDI(xp,F,InitParVec1,InitParVec2,gridSpacing,x,conn,nCelly)

    ParVec1 = F*InitParVec1;
    ParVec2 = F*InitParVec2;

    xEx = zeros(2,4);
    xEx(:,1) = [xp(1)-ParVec1(1)-ParVec2(1);xp(2)-ParVec1(2)-ParVec2(2)];
    xEx(:,2) = [xp(1)+ParVec1(1)-ParVec2(1);xp(2)+ParVec1(2)-ParVec2(2)];
    xEx(:,3) = [xp(1)+ParVec1(1)+ParVec2(1);xp(2)+ParVec1(2)+ParVec2(2)];
    xEx(:,4) = [xp(1)-ParVec1(1)+ParVec2(1);xp(2)-ParVec1(2)+ParVec2(2)];

    ie1 = (ceil(xEx(1,1)/gridSpacing))*nCelly-floor(xEx(2,1)/gridSpacing);
    ie2 = (ceil(xEx(1,2)/gridSpacing))*nCelly-floor(xEx(2,2)/gridSpacing);
    ie3 = (ceil(xEx(1,3)/gridSpacing))*nCelly-floor(xEx(2,3)/gridSpacing);
    ie4 = (ceil(xEx(1,4)/gridSpacing))*nCelly-floor(xEx(2,4)/gridSpacing);

    gridNodes = [conn(:,ie1);conn(:,ie2);conn(:,ie3);conn(:,ie4)];

    xpe1 = [2*(xEx(1,1)-x(1,conn(1,ie1)))/gridSpacing-1;2*(xEx(2,1)-x(2,conn(1,ie1)))/gridSpacing-1];
    xpe2 = [2*(xEx(1,2)-x(1,conn(1,ie2)))/gridSpacing-1;2*(xEx(2,2)-x(2,conn(1,ie2)))/gridSpacing-1];
    xpe3 = [2*(xEx(1,3)-x(1,conn(1,ie3)))/gridSpacing-1;2*(xEx(2,3)-x(2,conn(1,ie3)))/gridSpacing-1];
    xpe4 = [2*(xEx(1,4)-x(1,conn(1,ie4)))/gridSpacing-1;2*(xEx(2,4)-x(2,conn(1,ie4)))/gridSpacing-1];

    phis1 = 1/4*[(1-xpe1(1))*(1-xpe1(2)) (1+xpe1(1))*(1-xpe1(2)) ...
                 (1+xpe1(1))*(1+xpe1(2)) (1-xpe1(1))*(1+xpe1(2))];
    phis2 = 1/4*[(1-xpe2(1))*(1-xpe2(2)) (1+xpe2(1))*(1-xpe2(2)) ...
                 (1+xpe2(1))*(1+xpe2(2)) (1-xpe2(1))*(1+xpe2(2))];
    phis3 = 1/4*[(1-xpe3(1))*(1-xpe3(2)) (1+xpe3(1))*(1-xpe3(2)) ...
                 (1+xpe3(1))*(1+xpe3(2)) (1-xpe3(1))*(1+xpe3(2))];
    phis4 = 1/4*[(1-xpe4(1))*(1-xpe4(2)) (1+xpe4(1))*(1-xpe4(2)) ...
                 (1+xpe4(1))*(1+xpe4(2)) (1-xpe4(1))*(1+xpe4(2))];

    vol1 = 4.*(ParVec1(1)*ParVec2(2)-ParVec1(2)*ParVec2(1));
    gphi = [phis1*(ParVec1(2)-ParVec2(2)) phis2*(ParVec1(2)+ParVec2(2)) ...
           -phis3*(ParVec1(2)-ParVec2(2)) -phis4*(ParVec1(2)+ParVec2(2))
            phis1*(ParVec2(1)-ParVec1(1)) phis2*(-ParVec1(1)-ParVec2(1)) ...
           -phis3*(ParVec2(1)-ParVec1(1)) -phis4*(-ParVec1(1)-ParVec2(1))]/vol1;

end


function ids = bspline_ctrlP_ids(ctrlPx,ctrlPy,NctrlPx,deg)

deg = deg+1;
ids = zeros(deg*deg,1);
idy = (ctrlPy'-1)*NctrlPx;

for i = 1:deg

    ids((i-1)*deg+1:i*deg)=ctrlPx+idy(i);

end

end


function sfs = bspline_ctrlP_sfs(N_vecx,N_vecy,deg)

deg = deg+1;
sfs = zeros(deg*deg,1);

for i = 1:deg

    sfs((i-1)*deg+1:i*deg)=N_vecx'*N_vecy(i);

end

end


function dsfs = bspline_ctrlP_dsfs(N_vecx,N_vecy,B_vecx,B_vecy,deg)

deg = deg+1;
dsfs = zeros(2,deg*deg);

for i=1:deg

    dsfs(:,(i-1)*deg+1:i*deg)=[B_vecx*N_vecy(i);N_vecx*B_vecy(i)];

end

end
    

function [phi,nodes] = bsCorner(x,deg,knotx,knoty,NctrlPx)

[N_vecx,B_vecx,ctrlPx] = Bspline(knotx,deg,x(1));
[N_vecy,B_vecy,ctrlPy] = Bspline(knoty,deg,x(2));
nodes = bspline_ctrlP_ids(ctrlPx,ctrlPy,NctrlPx,deg);
phi = bspline_ctrlP_sfs(N_vecx,N_vecy,deg)';

end


function [x1,x2,x3,x4,ParVec1,ParVec2] = bsParDom(ParVec1,ParVec2,knotx,knoty,xp)

global gridSpacing

[x1,x2,x3,x4] = ParDomCorners(xp,ParVec1,ParVec2);
dumNx = size(knotx,2);
dumNy = size(knoty,2);
kntx1 = knotx(1)+1e-10*gridSpacing;
kntx2 = knotx(dumNx)-1e-10*gridSpacing;
knty1 = knoty(1)+1e-10*gridSpacing;
knty2 = knoty(dumNy)-1e-10*gridSpacing;

if min([x1(1) x2(1) x3(1) x4(1)])>=kntx1 && max([x1(1) x2(1) x3(1) x4(1)])<=kntx2 && min([x1(2) x2(2) x3(2) x4(2)])>=knty1 && max([x1(2) x2(2) x3(2) x4(2)])<=knty2

    return;

end

if xp(1)<kntx1 || xp(1)>kntx2 || xp(2)<knty1 || xp(2)>knty2
    
    disp('**** ERRORRRR - Particle outside the background domain');

end

if x1(1)<kntx1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(1)-ParVec2(1),1,kntx1,xp,ParVec1,ParVec2);

end

if x2(1)<kntx1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(1)-ParVec2(1),1,kntx1,xp,ParVec1,ParVec2);

end

if x3(1)<kntx1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(1)+ParVec2(1),1,kntx1,xp,ParVec1,ParVec2);

end

if x4(1)<kntx1

    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(1)+ParVec2(1),1,kntx1,xp,ParVec1,ParVec2);

end

if x1(1)>kntx2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(1)-ParVec2(1),1,kntx2,xp,ParVec1,ParVec2);

end

if x2(1)>kntx2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(1)-ParVec2(1),1,kntx2,xp,ParVec1,ParVec2);

end

if x3(1)>kntx2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(1)+ParVec2(1),1,kntx2,xp,ParVec1,ParVec2);

end

if x4(1)>kntx2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(1)+ParVec2(1),1,kntx2,xp,ParVec1,ParVec2);

end

if x1(2)<knty1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(2)-ParVec2(2),2,knty1,xp,ParVec1,ParVec2);

end

if x2(2)<knty1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(2)-ParVec2(2),2,knty1,xp,ParVec1,ParVec2);

end

if x3(2)<knty1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(2)+ParVec2(2),2,knty1,xp,ParVec1,ParVec2);

end

if x4(2)<knty1
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(2)+ParVec2(2),2,knty1,xp,ParVec1,ParVec2);

end

if x1(2)>knty2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(2)-ParVec2(2),2,knty2,xp,ParVec1,ParVec2);

end

if x2(2)>knty2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(2)-ParVec2(2),2,knty2,xp,ParVec1,ParVec2);

end

if x3(2)>knty2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(+ParVec1(2)+ParVec2(2),2,knty2,xp,ParVec1,ParVec2);

end

if x4(2)>knty2
    
    [x1,x2,x3,x4,ParVec1,ParVec2]=bsParDomSub(-ParVec1(2)+ParVec2(2),2,knty2,xp,ParVec1,ParVec2);

end

end


function [x1,x2,x3,x4,ParVec1,ParVec2] = bsParDomSub(r,i,knt,xp,ParVec1,ParVec2)

global gridSpacing

if abs(r)<1e-10*gridSpacing
    
    disp('**** ERRORRRR - Zero particle domain');

end

alpha = (knt-xp(i))/r;
ParVec1 = alpha*ParVec1;
ParVec2 = alpha*ParVec2;
[x1,x2,x3,x4] = ParDomCorners(xp,ParVec1,ParVec2);

end


function [xEx1,xEx2,xEx3,xEx4] = ParDomCorners(xp,ParVec1,ParVec2)

xEx1 = [xp(1)-ParVec1(1)-ParVec2(1);xp(2)-ParVec1(2)-ParVec2(2)];
xEx2 = [xp(1)+ParVec1(1)-ParVec2(1);xp(2)+ParVec1(2)-ParVec2(2)];
xEx3 = [xp(1)+ParVec1(1)+ParVec2(1);xp(2)+ParVec1(2)+ParVec2(2)];
xEx4 = [xp(1)-ParVec1(1)+ParVec2(1);xp(2)-ParVec1(2)+ParVec2(2)];
 
end