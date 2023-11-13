
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MPM-MATLAB code for CPDI, BSMPM, BSCPDI simulations                   %
%  Example 5.2 in 'Sadeghirad, B-spline convected particle domain       %
%   interpolation method, Engineering Analysis with Boundary Elements'  %
%                                                                       %
% Alireza Sadeghirad, PhD                                               %
%  sadeghirad@gmail.com                                                 %
%  November 2023                                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear
close all
global interpolator

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%       Preprocess        %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
interpolator = 'BSCPDI'; % input 'CPDI', 'BSMPM', or 'BSCPDI'
deg = 2;                 % degree of B-spline basis functions 
refinement_ratio = 1;    % input 1, 2, 4, 8, 16, ... to perform convergence study
Preprocess

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%  Loop over time steps  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for n = 1:totTimeSteps
    
    % update particle density & volume
    for i = 1:numpar
        detF = det(F(:,:,i));
        vol(i) = detF*vol0(i);
        rhop(i) = 1/detF*rho0(i);
    end

    % particle to grid: mass
    m = zeros(1,numnod);
    for i = 1:numpar
        for j = 1:nn(i)
            m(1,nodes(j,i)) = m(1,nodes(j,i))+phi(i,j)*rhop(i)*vol(i);
        end
    end

    % particle to grid: momentum
    q = zeros(2,numnod);
    for i = 1:numpar
        for j = 1:nn(i)
            q(:,nodes(j,i)) = q(:,nodes(j,i))+phi(i,j)*Vp(:,i)*rhop(i)*vol(i);
        end
    end

    % calculate external forces
    fext = zeros(2,numnod);
    for i = 1:numpar
        b_f = AnalyticalSolution(Xp(1,i),Xp(2,i),(n-1)*dt,'BodyForce');
        for j = 1:nn(i)
            fext(:,nodes(j,i)) = fext(:,nodes(j,i))+phi(i,j)*b_f*rhop(i)*vol(i);
        end
    end

    % calculate internal forces
    fint = zeros(2,numnod);
    for i = 1:numpar
        for j = 1:nn(i)
            fint(:,nodes(j,i)) = fint(:,nodes(j,i))+([gphi(1,j,i)*sigma(1,i)+gphi(2,j,i)*sigma(3,i);gphi(2,j,i)*sigma(2,i)+gphi(1,j,i)*sigma(4,i)])*vol(i);
        end
    end

    % solve equations of motion
    q1 = fext-fint;

    % enforce displacement BCs
    q(1,dbcx) = 0;
    q1(1,dbcx) = 0;
    q(2,dbcy) = 0;
    q1(2,dbcy) = 0;
    for i = dbcx_symm
        if abs(m(i(2)))>1e-12
            q(1,i(1)) = q(1,i(2))/m(i(2))*m(i(1));
            q1(1,i(1)) = q1(1,i(2))/m(i(2))*m(i(1));
        else
            q(1,i(1)) = 0;
            q1(1,i(1)) = 0;
        end
    end
    for i = dbcy_symm
        if abs(m(i(2)))>1e-12
            q(2,i(1)) = q(2,i(2))/m(i(2))*m(i(1));
            q1(2,i(1)) = q1(2,i(2))/m(i(2))*m(i(1));
        else
            q(2,i(1)) = 0;
            q1(2,i(1)) = 0;
        end
    end

    % update momentum at grid nodes
    q = q+q1*dt;
    
    % update velocity & position of particles
    dVp = zeros(2,numpar);
    dxp = zeros(2,numpar);
    for i = 1:numpar
        count = 1;
        while count<=nn(i)
            i5 = nodes(count,i);
            if abs(m(i5))>1e-10*mass_of_single_cell
                dVp(:,i) = dVp(:,i)+phi(i,count)*(q1(:,i5)/m(i5))*dt;
                dxp(:,i) = dxp(:,i)+phi(i,count)*(q(:,i5)/m(i5))*dt;
            end
            count = count+1;
        end
    end
    Vp = Vp+dVp;
    xp = xp+dxp;

    % update stress at particles
    sigma = zeros(4,numpar);
    for i = 1:numpar
        count = 1;
        vgn = zeros(2,nng(i));
        while count<=nng(i)
            i5 = gnodes(count,i);
            if abs(m(i5))>1e-10*mass_of_single_cell
                vgn(:,count) = q(:,i5)/m(i5);
            else
                vgn(:,count) = [0;0];
            end
            count = count+1;
        end
        LL = [gphi(1,1:nng(i),i)*vgn(1,:)' gphi(2,1:nng(i),i)*vgn(1,:)'
              gphi(1,1:nng(i),i)*vgn(2,:)' gphi(2,1:nng(i),i)*vgn(2,:)'];
        incrF = eye(2)+LL*dt;
        F(:,:,i) = incrF*F(:,:,i);
        Ssigma = Material(F(:,:,i),'Stress');
        sigma(:,i) = [Ssigma(1,1);Ssigma(2,2);Ssigma(1,2);Ssigma(2,1)];
    end

    % calculate effecte basis functions and their gradients
    [nn,nodes,phi,nng,gnodes,gphi] = CalcSFGradSF(x,nCelly,conn,xp,F,InitParVec1,InitParVec2,deg,knotx,knoty,NctrlPx,snode);

    % calculate error norms
    unum = xp-Xp;
    numpar = size(Vp,2);
    for i = 1:numpar
        uex = AnalyticalSolution(Xp(1,i),Xp(2,i),n*dt,'Displacement');
        u = unum(:,i)-uex;
        SumErrorDisp = SumErrorDisp+(u(1)*u(1)+u(2)*u(2))*vol(i);
        SumErrorDispEx = SumErrorDispEx+(uex(1)*uex(1)+uex(2)*uex(2))*vol(i);
        FzEx = AnalyticalSolution(Xp(1,i),Xp(2,i),(n)*dt,'DeformationGrad');
        Fz = F(:,:,i)-FzEx+eye(2);
        pe_ex = Material(FzEx,'Potential');
        pe_num = Material(Fz,'Potential');
        SumErrorEner = SumErrorEner+pe_num*vol(i);
        SumErrorEnerEx = SumErrorEnerEx+pe_ex*vol(i);
    end

    % plot snapshots
    figure(1)
    if n==round(totTimeSteps/4) || n==round(totTimeSteps/2) || n==round(totTimeSteps*3/4) || n==totTimeSteps
        subplot(2,2,round(n/totTimeSteps*4))
        Snapshot
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%   Print eroor norms  %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ee = (SumErrorEner/SumErrorEnerEx)^0.5;
Eu = (SumErrorDisp/SumErrorDispEx)^0.5;
disp('Displacement error norm:')
format long
disp(Eu)
disp('Energy error norm:')
disp(Ee)
