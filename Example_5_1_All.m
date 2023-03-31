clc; 
clear; 
close all;

%% Example 5.1 :
% Versteeg, H.K., Malalasekera, W., 2007. An introduction to computational 
% fuid dynamics: the finite volume method. Pearson Education. pp. 147-148

%% Notes:
% The CentralDifferencing, Upwind and QUICK differencing scheme have been 
% used to discretized the equations while the Gauss-Siedel iteration 
% method to solve the the set of algebraic equations.

%% Inputs

N=5;           % Number of nodes
ConvCrit=1e-4;  % Convergence criteria (for the Gauss-Seidel Scheme)
L=1.0;          % Length [m]
dx=L/N;         % Grid size [m]
rho=1.0;        % Density [kg m^-3]
u=2.5;          % Velocity [m s^-1]
F=rho*u;        % Convective flux term [kg m^-2 s^-1]
Gamma=0.1;      % Diffusion coefficient [kg m^-1 s^-1]
D=Gamma/dx;     % Diffusion conductance at cell faces [kg m^-2 s^-1]
Pe=F/D;         % Peclet number

%% Analytical solution %%

[distance_ana, phi_exact] = analytical_solution(L,rho,u,Gamma);

%%  Numerical solver using the FVM %%

[A_Upwind, b_Upwind] = discretisation(N,D,F,'Upwind');

[A_QUICK, b_QUICK] = discretisation(N,D,F,'QUICK');

[A_CentralDifferencing, b_CentralDifferencing] = ...
    discretisation(N,D,F,'CentralDifferencing');

%% Numerical Solution Using the FVM  %%

[x_Upwind, residual_Upwind, numItr_Upwind, phi_Upwind] = ...
    solution(N, ConvCrit, A_Upwind, b_Upwind);

[x_CentralDifferencing, residual_CentralDifferencing, ...
    numItr_CentralDifferencing, phi_CentralDifferencing] = ...
    solution(N, ConvCrit, A_CentralDifferencing, b_CentralDifferencing);

[x_QUICK, residual_QUICK, numItr_QUICK, phi_QUICK] = ...
    solution(N, ConvCrit, A_QUICK, b_QUICK);

distance_num=[dx/2:dx:L-dx/2];

%% Numerical solution deviations

RMSD_CentralDifferencing = RMSD(distance_num,phi_CentralDifferencing,rho,u,L,Gamma,N);
RMSD_Upwind = RMSD(distance_num,phi_Upwind,rho,u,L,Gamma,N);
RMSD_QUICK = RMSD(distance_num,phi_QUICK,rho,u,L,Gamma,N);

MaxD_CentralDifferencing = MaxD(distance_num,phi_CentralDifferencing,rho,u,L,Gamma);
MaxD_Upwind = MaxD(distance_num,phi_Upwind,rho,u,L,Gamma);
MaxD_QUICK = MaxD(distance_num,phi_QUICK,rho,u,L,Gamma);

%% Plot data

figure;

hold on
grid on

plot(distance_ana, phi_exact,'-k','LineWidth',1.5)
plot(distance_num,phi_CentralDifferencing,':sqr','LineWidth',1.5,'MarkerFaceColor','r')
plot(distance_num,phi_Upwind,':sqb','LineWidth',1.5,'MarkerFaceColor','b')
plot(distance_num,phi_QUICK,':sqg','LineWidth',1.5,'MarkerFaceColor','g')

set(gcf,'Units','centimeters');
afFigurePosition = [15 10 10 7.5];
set(gcf, 'Position', afFigurePosition);
set(gca,'xlim',[0 1],'xtick',[0:0.2:1.0],'FontSize',8,'FontWeight','normal');
% set(gca,'ylim',[0 1.1],'ytick',[0:0.1:1.1],'FontSize',8,'FontWeight','normal');
set(gcf,'color','w');
xlabel('Distance (m)','Fontsize',10); 
ylabel('$\phi$','interpreter','latex','FontSize',10);

lgd = legend('Exact solution','Numerical Central Differencing', ...
    'Numerical Upwind','Numerical QUICK');
lgd.FontSize = 10;

title(['Book example 5.1 with U = ',num2str(u),' m.s^{-1}'], ...
        'FontWeight','normal','fontsize',10);

dim = [.15 .3 .2 .1];
str = {'Inputs and numerical solutions deviations:','',...
    append('Number of nodes = ',num2str(N)), ...
    append('Velocity = ',num2str(u),' m.s^{-1}'), ...
    append('Peclet Number = ',num2str(Pe)), ...
    append('Central Differencing: RMSD = ', ...
    num2str(RMSD_CentralDifferencing),' MaxD = ', ...
    num2str(MaxD_CentralDifferencing)), ...
    append('Upwind: RMSD = ',num2str(RMSD_Upwind),' MaxD = ', ...
    num2str(MaxD_Upwind)), ...
    append('QUICK: RMSD = ',num2str(RMSD_QUICK),' MaxD = ', ...
    num2str(MaxD_QUICK))};
annotation('textbox',dim,'String',str,'FitBoxToText','on')

%% Functions

function [distance_ana, phi_exact] = analytical_solution(L,rho,u,Gamma)
    % Ua=1; Ub=0; % Boundary Conditions
    
    N2=100; % Number of nodes analytical solution
    
    distance_ana=zeros(N2+2,1);
    distance_ana(1,1)=0; 
    distance_ana(end,1)=L;
    
    phi_exact=zeros(N2+2,1);
    phi_exact(1,1)=1;
    phi_exact(end,1)=0;
    
    % Inner Points:
    
    Xa=L/N2;
    dxa=L/N2;
    
    AA=rho*u/Gamma;
    BB=rho*u*L/Gamma;
    
    for r=2:N2+1 % loop over a inner points
        phi_exact(r,1)=1+(1-exp(AA*Xa))/(exp(BB)-1);
        distance_ana(r,1)=Xa;
        Xa=Xa+dxa;
    end
end

function [A, b] = discretisation(N,D,F,scheme)
    switch scheme
        case 'CentralDifferencing'
            % Creating matrix A
    
            % Inner nodes:
            
            Sp=0;
            Su=0;
            ae=D-F/2;
            aw=D+F/2;
            ap=aw+ae;
            
            A=eye(N,N)*ap+diag(ones(1,N-1)*(-aw),-1)+diag(ones(1,N-1)*(-ae),1);
            
            % Assign boundary conditions at first and last node
            
            % First node:
            
            Sigma_A=1;
            
            Sp_A=-(2*D+F);
            Su_A=(2*D+F)*Sigma_A;
            
            aw = 0; 
            ae = D - F/2;
            
            ap=aw+ae-Sp_A;
    
            A(1,1) = ap;
            A(1,2) = -ae; 
            
            % Last node:
            
            Sigma_B=0;
    
            Sp_B=-(2*D-F);
            Su_B=(2*D-F)*Sigma_B;
            
            ae=0;
            aw=D+F/2;
            
            ap=aw+ae-Sp_B;
            
            A(N,N) = ap;
            A(N,N-1) = -aw;
            
            % Creating vector b:
            
            b=zeros(N,1);
            b(1,1)=Su_A;
            b(N,1)=Su_B;
    
        case 'Upwind'
            % Creating matrix A
    
            % Inner nodes:
            
            Sp=0;
            Su=0;
            ae=D+max(0,-F); % Note, Fw=Fe=F
            aw=D+max(F,0);
            ap=aw+ae-Sp;
            
            A=eye(N,N)*ap+diag(ones(1,N-1)*(-aw),-1)+diag(ones(1,N-1)*(-ae),1);
            
            % Assign boundary conditions at first and last node
            
            % First node:
            
            Sigma_A=1; % at x=0 (boundary condition)
            
            % source terms due to boundary conditions
            
            Sp=-(2*D+F); 
            Su_A=(2*D+F)*Sigma_A;
            
            % coefficients: 
            
            aw=0; 
            ap=aw+ae-Sp;
            A(1,1)=ap; % change in matrix A
            
            % Last node:
            
            Sigma_B=0;
            Sp=-(2*D);
            Su_B=(2*D)*Sigma_B;
            ae=0;
            aw=D+F;
            ap=aw+ae-Sp;
            A(N,N)=ap; % change in matrix A
            
            % Creating vector b:
            
            b=zeros(N,1);
            b(1,1)=Su_A;
            b(N,1)=Su_B;
    
        case 'QUICK'
            % Creating matrix A
    
            if F>0
                alpha = 1;
            else 
                alpha = 0;
            end 
            
            % Inner nodes (except for 1, 2 and N):
            
            Sp = 0;
            Su = 0;
            
            % Note, Fw=Fe=F
            
            ae = D - (3/8)*alpha*F - (6/8)*(1-alpha)*F - (1/8)*(1-alpha)*F; 
            aw = D + (6/8)*alpha*F + (1/8)*alpha*F + (3/8)*(1-alpha)*F;
            aww = -(1/8)*alpha*F;
            aee = (1/8)*(1-alpha)*F;
            
            ap = aw + ae + aww + aee - Sp;
            
            A = eye(N,N)*ap + diag(ones(1,N-1)*(-aw),-1) + diag(ones(1,N-1)*(-ae),1) + diag(ones(1,N-2)*(-aww),-2) + diag(ones(1,N-2)*(-aee),2);
            
            % Assign boundary conditions at first, second and last node
    
            % First node:
            
            Sigma_A = 1;
            
            Sp_A=-((8/3)*D + (2/8)*F + F);
            Su_A=((8/3)*D + (2/8)*F + F) * Sigma_A;
            
            aww = 0;
            aw = 0;
            ae = D + (1/3)*D - (3/8)*F;
            aee = 0;
            
            ap = aw + ae + aww - Sp_A;
            
            A(1,1) = ap;
            A(1,2) = -ae;
            A(1,3) = -aee;
            
            % Second node:
            
            Sp_B = (1/4)*F;
            Su_B = -(1/4)*F*Sigma_A;
            
            aww = 0;
            aw = D + (7/8)*F + (1/8)*F;
            ae = D - (3/8)*F;
            aee = 0;
            
            ap = aw + ae + aww - Sp_B;
            
            A(2,1) = -aw;
            A(2,2) = ap;
            A(2,3) = -ae;
            A(2,4) = -aee;
            
            % Last node:
            
            Sigma_C = 0;
            
            Sp_C = -((8/3)*D - F);
            Su_C = ((8/3)*D-F) * Sigma_C;
            
            aww = -(1/8)*F;
            aw = D + (1/3)*D + (6/8)*F;
            ae = 0;
            aee = 0;
            
            ap = aw + ae + aww - Sp_C;
            
            A(N,N-2) = -aww;
            A(N,N-1) = -aw;
            A(N,N) = ap;
            
            
            % Creating vector b:
            
            b=zeros(N,1);
            b(1,1)=Su_A; % Assign source term (such that Eq. 5.34 is correct)
            b(2,1)=Su_B;
            b(N,1)=Su_C;
            
            % Note that b(N,N)=Su_B=0
    
       otherwise
      
    end
end

function [x, residual, numItr, phi] = solution(N,ConvCrit,A,b)
    x0=zeros(N,1); % Initial guess of phi for the internal nodes
    
    % Gauss-Siedel Method for Solving Ax=b: 
    
    [x, residual, numItr] = gauss_seidel(A, b, x0, ConvCrit);
    
    phi=x; % The transported scalar 
end

function RMSD = RMSD(distance_num,phi,rho,u,L,Gamma,N)
    AA=rho*u/Gamma;
    BB=rho*u*L/Gamma;

    diff = (1+(1-exp(AA*(distance_num)))/(exp(BB)-1))' - phi;
    summation = sum(diff.^2);
    RMSD = ((1/N)*(summation))^(1/2);
end

function MaxD = MaxD(distance_num,phi,rho,u,L,Gamma)
    AA=rho*u/Gamma;
    BB=rho*u*L/Gamma;

    diff = (1+(1-exp(AA*(distance_num)))/(exp(BB)-1))' - phi;
    MaxD=max(abs(diff));
end
