
global gridSpacing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%       Input Data        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ncellx = 8;           % number of cells in x direction
ncelly = 8;           % number of cells in y direction
dt = 5e-4;            % time step
gridSpacing = 0.125;  % grid spacing
Ttot = 0.02;          % total simulation time
rho = 1000;           % mass density
ppc = 4;              % number of particles per cell

dt=dt/refinement_ratio;
gridSpacing=gridSpacing/refinement_ratio;
totTimeSteps=round(Ttot/dt);
mass_of_single_cell = rho*gridSpacing^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Determine size of SFs and gradSFs arrays %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(interpolator,'CPDI')

    snode = 16;

elseif strcmp(interpolator,'BSMPM') 

    snode = (deg+1)*(deg+1);

elseif strcmpi(interpolator,'BSCPDI')

    snode = 4*(deg+1)*(deg+1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     Generate the background grid     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCellx = ncellx*refinement_ratio+2;
nCelly = ncelly*refinement_ratio+2;

if strcmp(interpolator,'BSMPM') || strcmp(interpolator,'BSCPDI')

    nCellx=nCellx-2;
    nCelly=nCelly-2;

end

numnod = (nCellx+1)*(nCelly+1);
numelem = nCellx*nCelly;
x = zeros(2,numnod);
conn = zeros(4,numelem);

for i = 1:nCellx+1

    for j = 1:nCelly+1

        x(:,((nCelly+1)*(i-1)+j)) = [gridSpacing*(i-1);nCelly*gridSpacing-gridSpacing*(j-1)];

    end

end

for j = 1:nCellx

    for i = 1:nCelly

        elemn = (j-1)*nCelly+i;
        conn(4,elemn) = elemn+(j-1);
        conn(1,elemn) = conn(4,elemn)+1;
        conn(2,elemn) = conn(1,elemn)+nCelly+1;
        conn(3,elemn) = conn(2,elemn)-1;

    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%     Generate the parametric grid     %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NctrlPx=0; NctrlPy=0; knotx=[]; knoty=[];

if strcmp(interpolator,'BSMPM') || strcmp(interpolator,'BSCPDI')

    NctrlPx = nCellx+deg; NctrlPy=nCelly+deg;
    numnod = NctrlPx*NctrlPy;
    knotx = [ones(1,deg)*x(1,1) x(1,1:nCelly+1:nCellx*(nCelly+1)+1) ones(1,deg)*x(1,nCellx*(nCelly+1)+1)];
    knoty = [ones(1,deg)*x(2,nCelly+1) x(2,nCelly+1:-1:1) ones(1,deg)*x(2,1)];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%        Generate the particles        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numpar = 0;
p1 = zeros(1,sqrt(ppc));
p2 = zeros(2,ppc);
ii = 0;
xpp = zeros(2,ppc*numelem);
cond = zeros(1,ppc*numelem);

for i = 1:sqrt(ppc)
    
    for j = 1:sqrt(ppc)
        
        p1(j) = (2*j-1)/(sqrt(ppc)*2);
        p2(:,j+sqrt(ppc)*(i-1)) = [p1(j);p1(i)];
        
    end
    
end

for i = 1:numelem
    
    for j = 1:ppc
        
        ii = ii+1;
        xpp(:,ii) = x(:,conn(1,i))+p2(:,j)*gridSpacing;
        
        if strcmpi(interpolator,'CPDI')
            
            if (xpp(1,ii)-gridSpacing)>0 && (xpp(2,ii)-gridSpacing)>0 && (xpp(1,ii)-gridSpacing)<1 && (xpp(2,ii)-gridSpacing)<1
                
                cond(1,ii) = 1;
                
            end
            
        else
            
            if xpp(1,ii)>0 && xpp(2,ii)>0 && xpp(1,ii)<1 && xpp(2,ii)<1
                
                cond(1,ii) = 1;
                
            end
            
        end
        
        if cond(1,ii) == 1
            
            numpar=numpar+1;
            
        end
        
    end
    
end

xp = zeros(2,numpar);
ii=0;
jj=0;

for i = 1:numelem
    
    for j = 1:ppc
        
        ii = ii+1;
        
        if cond(1,ii) == 1
            
            jj = jj+1;
            xp(:,jj) = xpp(:,ii);
            
        end
        
    end
    
end

vol0 = gridSpacing*gridSpacing/ppc*ones(1,numpar);
InitParVec1 = [gridSpacing/(2*sqrt(ppc))*ones(1,numpar);zeros(1,numpar)];
InitParVec2 = [zeros(1,numpar);gridSpacing/(2*sqrt(ppc))*ones(1,numpar)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%         Displacement BCs        %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(interpolator,'BSMPM') || strcmp(interpolator,'BSCPDI')
    
    dbcx = [NctrlPx:NctrlPx:numnod 1:NctrlPx:numnod];
    dbcy = [(NctrlPy-1)*NctrlPx+1:numnod 1:NctrlPx];
    dbcx_symm = [];
    dbcy_symm = [];
    
elseif strcmp(interpolator,'CPDI')
    
    dbcy_u = find(abs(x(2,:))>gridSpacing*(nCelly-1)-1e-10*gridSpacing);
    dum1 = find(abs(x(2,:))>gridSpacing*nCelly-1e-10*gridSpacing);
    dum2 = setdiff(dbcy_u,dum1);
    dbcx_us = [dum1;dum2];
    dbcy_d = find(abs(x(2,:))<gridSpacing+1e-10*gridSpacing);
    dum1 = find(abs(x(2,:))<1e-10*gridSpacing);
    dum2 = setdiff(dbcy_d,dum1);
    dbcx_ds = [dum1;dum2];
    dbcx_r = find(abs(x(1,:))>gridSpacing*(nCellx-1)-1e-10*gridSpacing);
    dum1 = find(abs(x(1,:))>gridSpacing*nCellx-1e-10*gridSpacing);
    dum2 = setdiff(dbcx_r,dum1);
    dbcy_rs = [dum1;dum2];
    dbcx_l = find(abs(x(1,:))<gridSpacing+1e-10*gridSpacing);
    dum1 = find(abs(x(1,:))<1e-10*gridSpacing);
    dum2 = setdiff(dbcx_l,dum1);
    dbcy_ls = [dum1;dum2];
    dbcx = [dbcx_r dbcx_l];
    dbcy = [dbcy_u dbcy_d];
    dbcx_symm = [dbcx_us dbcx_ds];
    dbcy_symm = [dbcy_rs dbcy_ls];
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%       Initial conditions        %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
F = zeros(2,2,numpar); 
rho0 = rho*ones(1,numpar); 
sigma = zeros(4,numpar); 
q = zeros(2,numnod);

for i = 1:numpar

    F(1,1,i) = 1;
    F(2,2,i) = 1;

end

Xp = xp;
[nn,nodes,phi,nng,gnodes,gphi] = CalcSFGradSF(x,nCelly,conn,xp,F,InitParVec1,InitParVec2,deg,knotx,knoty,NctrlPx,snode);

if strcmp(interpolator,'CPDI')

    for i = 1:numnod

        q(:,i) = AnalyticalSolution(x(1,i),x(2,i),0,'Velocity');

    end

elseif strcmp(interpolator,'BSMPM') || strcmp(interpolator,'BSCPDI')

    numpar = size(Xp,2); 
    vp = zeros(2,numpar);
    actCtrlPs = zeros(1,numnod);
    
    for i = 1:size(nodes,2)

        for j = 1:nn(i)

            if phi(i,j)>1e-6

                actCtrlPs(nodes(j,i)) = 1; 

            end

        end

    end

    NactCtrlPs = sum(actCtrlPs);
    g2lCtrlP = zeros(1,numnod); 
    l2gCtrlP = zeros(1,NactCtrlPs); 
    count=0;

    for i = 1:numnod

        if actCtrlPs(i) == 1

            count = count+1; 
            g2lCtrlP(i) = count;
            l2gCtrlP(count) = i;

        end

    end

    A = sparse(numpar,NactCtrlPs); 

    for i = 1:numpar

        for j = 1:nn(i)

            if phi(i,j)>1e-6

                A(i,g2lCtrlP(nodes(j,i))) = A(i,g2lCtrlP(nodes(j,i)))+phi(i,j); 

            end

        end

    end

    for i = 1:numpar

        vp(:,i) = AnalyticalSolution(Xp(1,i),Xp(2,i),0,'Velocity');

    end

    VP1 = A'*vp(1,:)'; 
    VP2 = A'*vp(2,:)'; 
    A = A'*A;
    tol = 1e-8; 
    maxit = 1000; 
    qx = pcg(A,VP1,tol,maxit); 
    qy = pcg(A,VP2,tol,maxit);
    q(:,l2gCtrlP) = [qx';qy'];

end

m = zeros(1,numnod);

for i = 1:numpar

    for j = 1:nn(i)

        m(1,nodes(j,i)) = m(1,nodes(j,i))+phi(i,j)*rho0(i)*vol0(i);

    end

end

for i = 1:numnod

    q(:,i) = q(:,i)*m(1,i); 

end

dxp = zeros(2,numpar);

for i = 1:numpar

    count = 1;

    while count<=nn(i)

        i5 = nodes(count,i);

        if abs(m(i5))>1e-10*mass_of_single_cell

            dxp(:,i) = dxp(:,i)+phi(i,count)*(q(:,i5)/m(i5))*dt;

        end

        count = count+1;

    end

end

Vp = dxp/dt;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vol = vol0; 
rhop = rho0;
SumErrorDisp = 0; 
SumErrorDispEx = 0; 
SumErrorEner = 0; 
SumErrorEnerEx = 0;
