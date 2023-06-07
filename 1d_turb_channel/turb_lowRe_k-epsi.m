%---------------------------------------------------------------
% ---------------Low Reynolds Number k-epsilon--------------------------------
%----------------------------------------------------------------
% example of how to read DNS data.
%
% Channel data at Re_delta=7890, Re_tau=395, Re_theta=700.
% All quantites are normalized by u_tau and nu unless stated otherwise; delta denotes
% the half-channel width.  Data compiled from an unpublished work of Kim (1989).
%
% See also Kim, Moin & Moser (1987, JFM, vol 177) and Mansour, Kim & Moin (1988, JFM, vol 194).
%
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
epsi_tilde(1)=10^(-5);
for j=2:nj-1
   U(j)=1;
   k(j)=10^(-5);
   epsi_tilde(j)=10^(-5);
   vist(j)=10^(-5);
   f_mu(j)=1;
end
k(nj)=k(nj-1);
U(nj)=U(nj-1);
epsi_tilde(nj)=epsi_tilde(nj-1);
epsi_tilde(1)=1;
vist(nj)=vist(nj-1);
f_mu(nj)=f_mu(nj-1);
epsi = ones(1,nj)*1 ;
%% model constants , Eq. 3.37 ,(Versteegh & Malalasekera,1995)
C_mu=0.09 ;
sigma_k=1.00 ;
sigma_epsilon=1.30 ;
C_1epsilon=1.44 ;
C_2epsilon=1.92 ;
C_T = sqrt(2) ;
C_k = 0.5 ;
C_e = -2*C_k ;
%%
tao_w = 1 ;
% dudy(1) = 1/viscos(2) ;
kappa=0.41;
error=1;
count=0;
max_error=0.00001;

while error > max_error 

count = count+1;
 for j = 2 : nj - 1
% Compute the velocities at faces     
    Un(j) = fny(j)*U(j + 1) + (1 - fny(j))*U(j);
    Us(j) = fsy(j)*U(j - 1) + (1 - fsy(j))*U(j);
% Compute the velocity gradient du/dy
    dudy(j) = (Un(j) - Us(j)) / D_yn(j);

%%wall friction vel.
    tao_w = 1;
    ustar = sqrt(tao_w/rho);
% Compute the k at faces         
    kn(j) = fny(j)*k(j+1) + (1 - fny(j))*k(j);
    ks(j) = fsy(j)*k(j-1) + (1 - fsy(j))*k(j);
% Compute the velocity gradient dkdy
    dkdy(j) = (kn(j) - ks(j)) / D_yn(j);

% Compute the epsilon at faces         
    epsi_tilden(j) = fny(j)*epsi_tilde(j+1) + (1 - fny(j))*epsi_tilde(j);
    epsi_tildes(j) = fsy(j)*epsi_tilde(j-1) + (1 - fsy(j))*epsi_tilde(j);
% Compute the velocity gradient de~/dy
    depsi_tildedy(j) = (epsi_tilden(j) - epsi_tildes(j)) / D_yn(j);
    
% Compute the turbulent viscosity at faces         
    vistn(j) = fny(j)*vist(j+1) + (1 - fny(j))*vist(j);
    vists(j) = fsy(j)*vist(j-1) + (1 - fsy(j))*vist(j);
    
% Compute the production term         
    Pk(j) = vist(j)*dudy(j)^2;
    
% additional terms
    D(j) = 2* viscos(j)* k(j)/(ync(j)^2);
    epsi(j) = epsi_tilde(j) + D(j) ;
    epsi(nj)=epsi(nj-1);
    
    T_t(j) = max ( k(j)/epsi_tilde(j) , C_T*sqrt(viscos(j)/epsi(j))) ;
    R_y (j) = sqrt(k(j))*ync(j)/viscos(j) ;
    R_lamda(j) = ync(j)/(sqrt(viscos(j)*T_t(j))) ;
    f_mu(j) = 1 - exp(-0.01*R_lamda(j)-0.0068*R_lamda(j)^3) ;
%   shear stress
    u1u2(j)=vist(j)*dudy(j);
% Compute the k/e at faces         
    kdiven(j) = fny(j)*(k(j+1)/epsi(j+1)) + (1 - fny(j))*(k(j)/epsi(j));
    kdives(j) = fsy(j)*(k(j-1)/epsi(j-1)) + (1 - fsy(j))*(k(j)/epsi(j));
% Compute the velocity gradient de~/dy
    dkedy(j) = (kdiven(j) - kdives(j)) / D_yn(j);
    
%%
E_k(j) = C_k*vist(j)*min( dkedy(j)*depsi_tildedy(j),0) ;
E_e(j) = C_e*vist(j)/( T_t(j)^2)*(dkedy(j)*dkdy(j)) ;

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
    SU_k(j)  = (Pk(j))* D_yn(j);
    SP_k(j) = +( epsi(j)-E_k(j))*D_yn(j)/ k(j); 
    if j == nj - 1
        aN_k(j) = 0;
    end
    aP_k(j) = aN_k(j) + aS_k(j) + SP_k(j) ;

 %% Dissipation, epsi discretization
    nu_ne(j) = viscos(j) + vistn(j) / sigma_epsilon;
    nu_se(j) = viscos(j) + vists(j) / sigma_epsilon;
    aN_e(j)  = nu_ne(j) / Ypn_delta(j);
    aS_e(j)  = nu_se(j) / Ysp_delta(j);
    SU_e(j)  =(C_1epsilon*Pk(j)/T_t(j))*D_yn(j);
    SP_e(j)  = ((C_2epsilon*epsi_tilde(j)+D(j)*exp(-(R_y(j)/80)^2))/T_t(j)-E_e(j))* D_yn(j)/epsi_tilde(j); 
    if j == 2
        aS_e(j) = 0;
    end
    if j == nj - 1
        aN_e(j) = 0;
    end
    aP_e(j) = aN_e(j) + aS_e(j) + SP_e(j);
        
 end
 
 %% Boundary conditions Explicitly
    epsi_tilde(1) = epsi_tilde(2); %2*viscos*k(2)/y_node(2)^2;
    epsi (1) = epsi(2) ;
    U(nj) = U(nj-1);
    epsi_tilde(nj) = epsi_tilde(nj-1);
    k(nj) = k(nj-1);
 
%% Under relaxation factor , urf
%     urf=1;
%     U_old=U;
%     k_old=k;
%     epsi_old=epsi;
%     vist_old=vist;
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
       vist(j)= C_mu *f_mu(j)* k(j)*T_t(j);
       yplus(j)=ustar*ync(j)/viscos(2);
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
        
        epsi_tilde(j) = (aN_e(j)*epsi_tilde(j+1) + aS_e(j)*epsi_tilde(j-1) + SU_e(j)) / aP_e(j);
    end
% Convergence criterian
   for j=2:nj-1
      
% compute residuals R
%     R_u(j)      = sqrt((aN_u(j)*U(j + 1) + aS_u(j)*U(j - 1) + S_u_mod(j) - aP_u_mod(j)*U(j))^2);
%     R_k(j)      = sqrt((aN_k(j)*k(j + 1) + aS_k(j)*k(j - 1) + S_k_mod(j) - aP_k_mod(j)*k(j))^2);
%     R_epsi(j)   = sqrt((aN_e(j)*epsi(j + 1) + aS_e(j)*epsi(j - 1) + S_e_mod(j) - aP_e_mod(j)*epsi(j))^2);
    
    R_u(j)      = sqrt((aN_u(j)*U(j + 1) + aS_u(j)*U(j - 1) + SU_u(j) - aP_u(j)*U(j))^2);
    R_k(j)      = sqrt((aN_k(j)*k(j + 1) + aS_k(j)*k(j - 1) + SU_k(j) - aP_k(j)*k(j))^2);
    R_epsi(j)   = sqrt((aN_e(j)*epsi_tilde(j + 1) + aS_e(j)*epsi_tilde(j - 1) + SU_e(j) - aP_e(j)*epsi_tilde(j))^2);

% compute the flux F
    F(j)= U(j)^2*D_yn(j);
   end
   error_u(count) = sum(R_u) / sum(F);
   error_k(count) = sum(R_k) / sum(F);
   error_e(count) = sum(R_epsi) / sum(F);
   error(count)   = max([error_u(count),error_k(count),error_e(count)]);
%   if ((count == 30000))
%        break
%   end
%    if ((k < 0))
%        break
%   end
   
%   figure(1)
% plot(u_dns(1:nj-1,1),yplus','b',U(1,1:nj-1),yplus,'r')
% figure(2)
% plot(k_dns,yplus,'bo',k(1,1:nj-1),yplus,'r')
% figure(3)
% eps_dns=dns_data(:,2)*ustar^4/viscos(2);
% plot(eps_dns, yplus,'bo',epsi(1,1:nj-1),yplus,'m')
% ylim([0,30])
% figure(4)
% plot(yplus,vist+viscos)
% xlim([0,100])
% figure(5)
% counti(count)=count;
% plot(counti)
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

% Compare also with the different terms in the k-eq. 
% Read DNS data from file 'dns_data.dat'
%
% 6 coulumns as below:
%
%      y+         Diss        prod     vel_p_grad   Turb_diff   Visc_diff
%
% Please note that all terms are normalized by ustar^4/viscos
%
%shear stress

figure(5)
plot(u1u2,ync(2:end),'r'); hold on
plot(-uv_dns,y_dns,'bo')
legend('k-epsi','DNS')
ylabel('y-cordinat')
xlabel('reynolds stress')
print reynoldsstress.eps -deps

% kolmogorov scales
for i= 1:length(epsi)
    k_l(i)=((viscos(i)^3)/epsi(i))^0.25;
    k_t(i)=(viscos(i)/epsi(i))^0.5;
    k_v(i)=(viscos(i)*epsi(i))^0.25;
end
% Large turbulent scales
l_scale=(R_lamda.^0.75)/10;
t_scale=(R_lamda.^0.5)/10;
v_scale=(R_lamda.^0.25)/10;

figure(6)
plot(k_l,ync,'b'); hold on
plot(k_t,ync,'r')
plot(k_v,ync,'g')
plot(l_scale,ync(2:98),'b--')
plot(t_scale,ync(2:98),'r--')
plot(v_scale,ync(2:98),'g--')
legend('Length scale','Time scale','Velocity scale','Large length scale/10',...
    'Large time scale/10','Large velocity scale/10','location','best')
ylabel('y-cordinat')
xlabel('[m] / [sec] / [m/sec]')
xlim([0 0.9])
print scales.eps -depsc

% Budget
for i=1:length(epsi)-1
    production(i)= u1u2(i)*dudy(i)*viscos(i);
end

vel_p_grad=viscos.*gradient(U.*(1/rho),ync);
turb_diff=viscos./sigma_k.*vist.*gradient(k,ync);
visc_diff=viscos.*gradient(gradient(k,ync),ync);

figure(7)
% dns values 
plot(dns_data(:,2),y_dns,'r--');hold on% epsi
plot(dns_data(:,3),y_dns,'g--')
plot(dns_data(:,4),y_dns,'b--')
plot(dns_data(:,5),y_dns,'m--')n
plot(dns_data(:,6),y_dns,'k--')
% computed 
plot(epsi.*viscos,ync,'r');   
plot(production,ync(2:98),'g')    
%plot(vel_p_grad,ync,'b')    
plot(turb_diff,ync,'m')     
%plot(visc_diff,ync,'k')     
ylim([0 0.5])
legend('Dissipation DNS','Production DNS','pv fluctuation DNS',...
    'Turbulent diffusion DNS','Viscous diffusion DNS','Dissipation Calc.',...
    'Production Calc.','Turbulent diffusion calc.')
ylabel('y-coordinate')
print budget.eps -depsc