
% K3
%
close all
clear all
format long

% wall friction velocity
ustar=1;
% read u from DNS data base
load u_dns.dat

% read y from DNS data base
load y_dns.dat
y=y_dns;          %yc is a vector contains the faces coordinates (grid)
nj=length(y)+1; % nj: number of node points  
u_dns(nj)=u_dns(nj-1);
load u2_dns.dat
load v2_dns.dat
load w2_dns.dat
load dns_data.dat
k_dns=0.5*(u2_dns+v2_dns+w2_dns);
% eps_dns=dns_data(:,2)*ustar^4/viscos(2);

%% grid properties
ync(1)=y(1);
for j=2:nj-1
    ync(j)  = (y(j)+y(j-1))/2; 
end
ync(nj)=y(nj-1);

for  j=2:nj-1
    Ypn_delta(j)    = ync(j + 1) - ync(j);
    Ysp_delta(j)    = ync(j) - ync(j - 1);
    D_yn(j)          = y(j) - y(j - 1);    
    fny(j)          = 0.5*D_yn(j) / Ypn_delta(j);
    fsy(j)          = 0.5*D_yn(j) / Ysp_delta(j);
end

%%
% viscosity
viscos=ones(1,nj)*1/395;
viscos(1)=0;
%density
rho = 1 ;
% init velocity, vist, k & eps
U(1)=0;
k(1)=0;
epsi(1)=10^(-5);
for j=2:nj-1
   U(j)=1;
   k(j)=10^(-5);
   epsi(j)=10^(-5);
   vist(j)=10^(-5);
end
k(nj)=k(nj-1);
U(nj)=U(nj-1);
epsi(nj)=epsi(nj-1);
epsi(1)=1;
vist(nj)=vist(nj-1);

%% model constants , Eq. 3.37 ,(Versteegh & Malalasekera,1995)
C_mu=0.09 ;
sigma_k=1.00 ;
sigma_epsilon=1.30 ;
C_1epsilon=1.44 ;
C_2epsilon=1.92 ;

%%
tao_w = 1 ;
% dudy(1) = 1/viscos(2) ;
kappa=0.41;
error=1;
count=0;
max_error=0.001;

while error > max_error 

count = count+1;
 for j = 2 : nj - 1
% Compute the velocities at faces     
    Un(j) = fny(j)*U(j + 1) + (1 - fny(j))*U(j);
    Us(j) = fsy(j)*U(j - 1) + (1 - fsy(j))*U(j);
% Compute the velocity gradient du/dy
    dudy(j) = (Un(j) - Us(j)) / D_yn(j);
%     dudy(j) = (U(j) - U(j-1)) / D_yn(j);

%%wall friction vel.
    tao_w = 1;
    ustar = sqrt(tao_w/rho);
% Compute the k at faces         
    kn(j) = fny(j)*k(j+1) + (1 - fny(j))*k(j);
    ks(j) = fsy(j)*k(j-1) + (1 - fsy(j))*k(j);
    
% Compute the epsilon at faces         
    epsin(j) = fny(j)*epsi(j+1) + (1 - fny(j))*epsi(j);
    epsis(j) = fsy(j)*epsi(j-1) + (1 - fsy(j))*epsi(j);
    
% Compute the turbulent viscosity at faces         
    vistn(j) = fny(j)*vist(j+1) + (1 - fny(j))*vist(j);
    vists(j) = fsy(j)*vist(j-1) + (1 - fsy(j))*vist(j);
    
% Compute the production term         
    Pk(j) = vist(j)*dudy(j)^2;

%% Mean velocity, U discretization
    nu_n(j) = viscos(j) + vistn(j);
    nu_s(j) = viscos(j) + vists(j);
    aN_u(j) = nu_n(j) / Ypn_delta(j);
    aS_u(j) = nu_s(j) / Ysp_delta(j);
    SU_u(j) = tao_w * D_yn(j);

    if j == nj-1
    aN_u(j) = 0;
    end
    aP_u(j) = aN_u(j) + aS_u(j);
 %% Turbulent kinetic energy, k discretization
    nu_nk(j) = viscos(j) + vistn(j) / sigma_k;
    nu_sk(j) = viscos(j) + vists(j) / sigma_k;
    aN_k(j)  = nu_nk(j) / Ypn_delta(j);
    aS_k(j)  = nu_sk(j) / Ysp_delta(j);
    SU_k(j)  = Pk(j) * D_yn(j);
    SP_k(j) = - epsi(j)*D_yn(j)/ k(j); 
    if j == nj - 1
        aN_k(j) = 0;
    end
    aP_k(j) = aN_k(j) + aS_k(j) - SP_k(j) ;

 %% Dissipation, epsi discretization
    nu_ne(j) = viscos(j) + vistn(j) / sigma_epsilon;
    nu_se(j) = viscos(j) + vists(j) / sigma_epsilon;
    aN_e(j)  = nu_ne(j) / Ypn_delta(j);
    aS_e(j)  = nu_se(j) / Ysp_delta(j);
    SU_e(j)  = C_1epsilon*(epsi(j) / k(j))*Pk(j)*D_yn(j);
    SP_e(j)  = -C_2epsilon*(epsi(j) / k(j))* D_yn(j); 
    if j == 2
        aS_e(j) = 0;
    end
    if j == nj - 1
        aN_e(j) = 0;
    end
    aP_e(j) = aN_e(j) + aS_e(j) - SP_e(j);
        
 end
 
 %% Boundary conditions Explicitly
    epsi(1) = epsi(2); %2*viscos*k(2)/y_node(2)^2;
    U(nj) = U(nj-1);
    epsi(nj) = epsi(nj-1);
    k(nj) = k(nj-1);
 
%% Under relaxation factor , urf
%     urf=1;
%     U_old=U;
%     k_old=k;
%     epsi_old=epsi;
    vist_old=vist;
% 
%     aP_u_mod=aP_u/urf;
%     aP_k_mod=aP_k/urf;
%     aP_e_mod=aP_e/urf;
%     for j = 2 : nj - 1
%         S_u_mod(j)=SU_u(j) + aP_u_mod(j)*(1-urf)*U_old(j);
%         S_k_mod(j)=SU_k(j) + aP_k_mod(j)*(1-urf)*k_old(j);
%         S_e_mod(j)=SU_e(j) + aP_e_mod(j)*(1-urf)*epsi_old(j);
%     end
%     
   if count < 2000   % use mixing length model for turbulent viscosity if count >2000
      for j=2:nj-1
% compute turbulent viscosity         
         yplus(j)=ustar*ync(j)/viscos(2);
         damp=1-exp(-yplus(j)/26);
         ell=min(damp*kappa*ync(j),0.09);
%          vist(j)=urf*abs(dudy(j-1))*ell^2+(1-urf)*vist_old(j);
         vist(j)=abs(dudy(j))*ell^2 ;

      end
   else
       for j=2:nj-1
       vist(j)= C_mu * k(j)^2 / epsi(j);
       end
   end
   vist(nj)=vist(nj-1);
%
%% Gauss-Seidel Solver
%     for j = 2 : nj - 1        
%         % Calculate U, k and epsilon
%         U_temp(j) = (aN_u(j)*U(j + 1) + aS_u(j)*U(j - 1) + S_u_mod(j)) / aP_u_mod(j);
%         U(j) = urf*U_temp(j) + (1-urf)*U_old(j);
%         
%         k_temp(j) = (aN_k(j)*k(j + 1) + aS_k(j)*k(j - 1) + S_k_mod(j)) / aP_k_mod(j);
%         k(j) = urf*k_temp(j) + (1-urf)*k_old(j);
%         
%         epsi_temp(j) = (aN_e(j)*epsi(j + 1) + aS_e(j)*epsi(j - 1) + S_e_mod(j)) / aP_e_mod(j);
%         epsi(j) = urf*epsi_temp(j) + (1-urf)*epsi_old(j);
%     end
   for j = 2 : nj - 1        
        % Calculate U, k and epsilon
        U(j) = (aN_u(j)*U(j+1) + aS_u(j)*U(j-1) + SU_u(j)) / aP_u(j);
        
        k(j) = (aN_k(j)*k(j+1) + aS_k(j)*k(j-1) + SU_k(j)) / aP_k(j);
        
        epsi(j) = (aN_e(j)*epsi(j+1) + aS_e(j)*epsi(j-1) + SU_e(j)) / aP_e(j);
        gamma(j) = k(j)/epsi(j) ;
    end
% Convergence criterian
   for j=2:nj-1
% compute residuals R
%     R_u(j)      = sqrt((aN_u(j)*U(j + 1) + aS_u(j)*U(j - 1) + S_u_mod(j) - aP_u_mod(j)*U(j))^2);
%     R_k(j)      = sqrt((aN_k(j)*k(j + 1) + aS_k(j)*k(j - 1) + S_k_mod(j) - aP_k_mod(j)*k(j))^2);
%     R_epsi(j)   = sqrt((aN_e(j)*epsi(j + 1) + aS_e(j)*epsi(j - 1) + S_e_mod(j) - aP_e_mod(j)*epsi(j))^2);
    
    R_u(j)      = sqrt((aN_u(j)*U(j + 1) + aS_u(j)*U(j - 1) + SU_u(j) - aP_u(j)*U(j))^2);
    R_k(j)      = sqrt((aN_k(j)*k(j + 1) + aS_k(j)*k(j - 1) + SU_k(j) - aP_k(j)*k(j))^2);
    R_epsi(j)   = sqrt((aN_e(j)*epsi(j + 1) + aS_e(j)*epsi(j - 1) + SU_e(j) - aP_e(j)*epsi(j))^2);

% compute the flux F
    F(j)= U(j)^2*D_yn(j);
   end
   error_u(count) = sum(R_u) / sum(F);
   error_k(count) = sum(R_k) / sum(F);
   error_e(count) = sum(R_epsi) / sum(F);
   error(count)   = max([error_u(count),error_k(count),error_e(count)]);
%   if count == 10000
%        break
%   end
   
end  %while

  figure(1)
plot(u_dns(1:nj-1,1),yplus','b',U(1,1:nj-1),yplus,'r')
% plot(u_dns,ync,'b',U,ync,'r')
figure(2)
plot(k_dns,yplus,'bo',k(1,1:nj-1),yplus,'r')
figure(3)
eps_dns=dns_data(:,2)*ustar^4/viscos(2);
plot(eps_dns, yplus,'bo',epsi(1,1:nj-1),yplus,'m')
figure(4)
plot(yplus,vist(1:nj-1)+viscos(1:nj-1))
