%% GTD Library Generation 
% Using coordinate change (omega_e1=0), normalisation (|omega|=1) and
% diagonalisation (purely imaginary), generate library that can be
% interpolated and recreat diffusivity and average swimming. 
clear all
% par=parpool(28);
%% Define discretization parameters
settings.n_phi=32; % Has to be even for FFT. 2^N recommended
settings.n_theta=101; % Had better to be 1+(multiples of 5,4,3 or 2) for Newton-Cote
settings.tn_phi=settings.n_phi/2+1;
settings.tn_theta=(floor(settings.n_theta/2)+1);
settings.N_total=(settings.n_theta)*settings.n_phi;
settings.dtheta=(pi/(settings.n_theta-1));
settings.dphi=2*pi/(settings.n_phi);
% N_total=settings.n_phi*settings.n_theta;
settings.theta=[0:settings.dtheta:pi]';
settings.phi=[0:settings.dphi:(2*pi-settings.dphi)];
settings.kphi=[0:(settings.tn_phi-1)];

p1_field=sin(settings.theta)*cos(settings.phi);
p2_field=sin(settings.theta)*sin(settings.phi);
p3_field=cos(settings.theta)*ones(size(settings.phi));

%% Newton-Cote Integrand
    if mod(settings.n_theta,5)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[19 repmat([75 50 50 75 38],1,(settings.n_theta-6)/5)  75 50 50 75 19]*5/288*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,4)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[7 repmat([32 12 32 14],1,(settings.n_theta-5)/4)  32 12 32 7]*2/45*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,3)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([3 3 2],1,(settings.n_theta-4)/3) 3 3 1]*3/8*settings.dtheta*2*pi;
    elseif mod(settings.n_theta,2)==1
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0].*[1 repmat([4 2],1,(N-3)/2) 4 1]/3*settings.dtheta*2*pi;
    else
        settings.integrad=[0 sin(settings.theta(2:end-1)') 0]*settings.dtheta*pi*2; %% Trapezoid Rule
    end

settings.int=kron(settings.integrad,ones(1,settings.n_phi)/settings.n_phi);
%% Define input parameters
settings.beta=0.; %gyrotactic 

% G, the TRANSPOSE of velocity gradient (u=x.G) OR Jacobian of velocity 
G11=0; % du/dx
G12=0; % dv/dx
G13=1; % dw/dx
G21=0; % du/dy
G22=0; % dv/dy
G23=0; % dw/dy
G31=0; % du/dz
G32=0; % dv/dz
G33=0; % dw/dz

settings.omega1=G23-G32;
settings.omega2=G31-G13;
settings.omega3=G12-G21;

settings.E11=G11;
settings.E12=(G12+G21)/2;
settings.E13=(G13+G31)/2;
settings.E22=G22;
settings.E23=(G23+G32)/2;
settings.E33=G33;

settings.B=0.9; % Bretherton constant

%% Initialisation
%S_loop=[0:1:100];
S_loop=2;
% G12_loop=[-exp(6:-0.025:-12) 0 exp(-12:0.025:6)];

% Looping Mesh
N_loop=numel(S_loop);

% Initialise result array
     D_array=zeros(N_loop,9);
    bf_array=zeros(N_loop,3);
 pavg1_array=zeros(N_loop,1);
 pavg2_array=zeros(N_loop,1);
 pavg3_array=zeros(N_loop,1);

%% Loop for different S
for ii=1:N_loop
    %% Set up parameters in loop
    sloc=settings;
%     sloc.S=S_mesh(ii);
    sloc.omega1=sloc.omega1*S_loop(ii);
    sloc.omega2=sloc.omega2*S_loop(ii);
    sloc.omega3=sloc.omega3*S_loop(ii);

    
    sloc.E11=sloc.E11*S_loop(ii);
    sloc.E12=sloc.E12*S_loop(ii);
    sloc.E13=sloc.E13*S_loop(ii);
    sloc.E22=sloc.E22*S_loop(ii);
    sloc.E23=sloc.E23*S_loop(ii);
    sloc.E33=sloc.E33*S_loop(ii);

    
    %% Solving for f(e)
    LHSL=L_FD_LHS(sloc);
    LHS_F=L_FD_LHS_BC(LHSL,sloc);
    
    % Invert LHS_F with zeros on right
    [V,~]=eigs(LHS_F,1,0);
    f_sol=V/(kron(sloc.integrad,ones(1,sloc.n_phi)/sloc.n_phi)*V);

    % Alternative way to invert LHS_F
%     LHS_F(sloc.N_total/2,:)=kron(sloc.integrad,ones(1,sloc.n_phi)/sloc.n_phi);
%     RHS_F=zeros(sloc.n_theta*sloc.n_phi,1);
%     RHS_F(sloc.N_total/2,1)=1;
%     f_sol=LHS_F\RHS_F;

    f0=transpose(reshape(f_sol,sloc.n_phi,sloc.n_theta));

    %% Finding p_avg 
    % Define weight function prefactor
    weight = - p3_field;
    pavg1_array(ii)=sloc.integrad*mean(f0.*p1_field.*weight,2);
    pavg2_array(ii)=sloc.integrad*mean(f0.*p2_field.*weight,2);
    pavg3_array(ii)=sloc.integrad*mean(f0.*p3_field.*weight,2);    
    
    %% b_c
    RHS_b1=f0.*(weight.*p1_field-pavg1_array(ii));RHS_b1([1 sloc.n_theta],:)=0;
    RHS_b2=f0.*(weight.*p2_field-pavg2_array(ii));RHS_b2([1 sloc.n_theta],:)=0;
    RHS_b3=f0.*(weight.*p3_field-pavg3_array(ii));RHS_b3([1 sloc.n_theta],:)=0;
    
    LHS_b=[LHS_F sloc.int';sloc.int 0];

    RHS_b1_col=reshape(transpose(RHS_b1),(sloc.n_theta)*sloc.n_phi,1);
    RHS_b2_col=reshape(transpose(RHS_b2),(sloc.n_theta)*sloc.n_phi,1);
    RHS_b3_col=reshape(transpose(RHS_b3),(sloc.n_theta)*sloc.n_phi,1);
    
    b1col=(-LHS_b)\[RHS_b1_col;0];
    b2col=(-LHS_b)\[RHS_b2_col;0];
    b3col=(-LHS_b)\[RHS_b3_col;0];
    
    %% f_v_s
    % Invert LHS_F with zeros on the right
    % Impose integral restriction
    LHS_f=[LHS_F sloc.int';sloc.int 0];
    RHS_f=zeros(sloc.n_theta*sloc.n_phi+1,1);
    [W,~]=eigs(LHS_f,1,0);
    f_vsol=W/([kron(sloc.integrad,ones(1,sloc.n_phi)/sloc.n_phi) 0]*W);
    f_v_s=transpose(reshape(f_vsol(1:3232),sloc.n_phi,sloc.n_theta));
    
    %% D_c
    D11_field=transpose(reshape(b1col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p1_field;
    D12_field=transpose(reshape(b2col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p1_field;
    D13_field=transpose(reshape(b3col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p1_field;
    D21_field=transpose(reshape(b1col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p2_field;
    D22_field=transpose(reshape(b2col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p2_field;
    D23_field=transpose(reshape(b3col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p2_field;
    D31_field=transpose(reshape(b1col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p3_field;
    D32_field=transpose(reshape(b2col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p3_field;
    D33_field=transpose(reshape(b3col(1:sloc.N_total),sloc.n_phi,sloc.n_theta)).*weight.*p3_field;
    
    D_temp=zeros(9,1);
    D_temp(1)=sloc.integrad*mean(D11_field,2);
    D_temp(2)=sloc.integrad*mean(D12_field,2);
    D_temp(3)=sloc.integrad*mean(D13_field,2);
    D_temp(4)=sloc.integrad*mean(D21_field,2);
    D_temp(5)=sloc.integrad*mean(D22_field,2);
    D_temp(6)=sloc.integrad*mean(D23_field,2);
    D_temp(7)=sloc.integrad*mean(D31_field,2);
    D_temp(8)=sloc.integrad*mean(D32_field,2);
    D_temp(9)=sloc.integrad*mean(D33_field,2);
    
    D_array(ii,:)=D_temp;
    
    
    disp([num2str(ii) '/' num2str(N_loop)]);

    
end
%
% figure(7)
% plot(S_loop,pavg1_array)
% grid on
% xlabel('S'), ylabel('(K_g)_1')
% print('-dpdf')
% %
% figure(8)
% plot(S_loop,pavg2_array)
% grid on
% xlabel('S'), ylabel('(K_g)_2')
% print('-dpdf')
% %
% figure(9)
% plot(S_loop,pavg3_array)
% grid on
% xlabel('S'), ylabel('(K_g)_3')
% print('-dpdf')

% res_array=reshape(res_array,numel(G12_loop),numel(G21_loop),numel(G11_loop),numel(G22_loop),NBvar); %#ok<*NASGU>
% e1_array=reshape(e1_array,numel(G12_loop),numel(G21_loop),numel(G11_loop),numel(G22_loop));
% e2_array=reshape(e2_array,numel(G12_loop),numel(G21_loop),numel(G11_loop),numel(G22_loop));
% e3_array=reshape(e3_array,numel(G12_loop),numel(G21_loop),numel(G11_loop),numel(G22_loop));

%% Saving
% name=['GTD_beta_' num2str(settings.beta*10) '_GTD_lib_' num2str(settings.n_phi) '_' num2str(settings.n_theta) '_Sbase_ELComp'];
% clearvars settings e_all_field e1_field e2_field e3_field par;
% save([name '.mat']);
% exit