%------------------------------------------------------------------------------------------
%                      2D FINITE VOLUME CODE FOR CONVECTION DIFFUSION (@Surya Kiran Peravali)
%                      CHALMERS UNIVERSITY OF TECHNOLOGY, Gothenburg,2014
%------------------------------------------------------------------------------------------

clc
clear all
close all
tic

T_A = 20 ; %inlet temp.
T_B = 10;
T_i = 0 ; %tempretature initialization
% solver flag
% for TDMA         ... 1
% for Gauss-Seidel ... 0
s_flag=0 ;
% boundary condition flag
% for Dirichlet at right wall ... 1
% for Neumann at right wall   ... 0
bc_flag = 1 ;
inlet_node = 4 ;
outlet_node = 4 ;
kfactor=10; %for increasing or decreasing k value
error = 0.001 ; % number of iterations
kn = 1 ;
c_p = 500 ;
gamma = kn*kfactor/c_p ;
rho = 1;
s_bar = 0; % source term

% read xc
load xc.dat
nim1=length(xc);
% nim1 = ni-1 = number of grid lines. Number of cell nodes = ni
ni=nim1+1;
%
% read yc
load yc.dat
njm1=length(yc);
nj=njm1+1;
%
% read u
load u.dat
u2d=reshape(u,ni,nj);
% read v
load v.dat
v2d=reshape(v,ni,nj);
%
% xc and yc are the coordinates of the grid
%
%       o---------------o  yc(j)
%       |               |
%       |               |
%       |       P       |
%       |               |
%       |               |
%       o---------------o  yc(j-1)
%
%     xc(i-1)          xc(i)          
%
% The cell-node (P) has the coordinates xp(i), yp(j)
%
% compute the x-coordinates of the cell centres
x = xc ;
y = yc ;

for i=2:nim1
   xnc(i)=0.5*(x(i)+x(i-1));
end
xnc(1)=x(1);
xnc(ni)=x(nim1);
%
% take the transpose of x
xnc=xnc';

% compute the y-coordinates of the cell centres
for j=2:njm1
   ync(j)=0.5*(y(j)+y(j-1));
end
ync(1)=y(1);
ync(nj)=y(njm1);
%
% take the transpose of y
ync=ync';

for i = 1:(length(x)-1)
    xn(i) = (x(i+1)- x(i))/2; %delta x/2
end
for j = 1:(length(y)-1)
    yn(j) = (y(j+1)-y(j))/2; %delta y/2
end

    for  k = 1:(ni-2)
Xpe_delta(k)= xnc(k+2)-xnc(k+1);
Xwp_delta(ni-k-1)= xnc(ni-k)-xnc(ni-k-1);
    end
    
    for  l = 1:(nj-2)
Ypn_delta(l)= ync(l+2)-ync(l+1);
Ysp_delta(nj-l-1)= ync(nj-l)-ync(nj-l-1);
    end
    
fex= xn./Xpe_delta;
fwx= xn./Xwp_delta;
fny = yn./Ypn_delta;
fsy = yn./Ysp_delta;
%%

u_f_x = zeros(nj-2,ni-1) ;
for s= 1:nj-2
        for z=1:ni-3
 u_f_x(s,z+1)=interp1(xnc,u2d(s+1,:),x(z+1))     ;
%  Tf_x(s,z+1)= (T(s+1,z+2)-T(s+1,z+1))*(x(z+1)-xnc(z))/(xnc(z+1)-xnc(z))+T(s+1,z+1) ;
        end
end
u_f_x(:,1)= u2d(2:nj-1,1);
u_f_x(:,ni-1)= u2d(2:ni-1,ni);
u_w = u_f_x(:,1:ni-2);
u_e = u_f_x(:,2:ni-1);

u_f_y = zeros(nj-1,ni-2) ;
for z= 1:ni-2
        for s=1:nj-3
 u_f_y(s+1,z)=interp1(ync,v2d(:,z+1),y(s+1)) ;
%  Tf_x(s,z+1)= (T(s+1,z+2)-T(s+1,z+1))*(x(z+1)-xnc(z))/(xnc(z+1)-xnc(z))+T(s+1,z+1) ;
        end
end
u_f_y(1,:)= v2d(1,2:ni-1);
u_f_y(nj-1,:)= v2d(nj,ni-1);

v_s = u_f_y(1:nj-2,:);
v_n = u_f_y(2:nj-1,:);
%%
gamma_node = ones(nj,ni)*gamma;
gamma_e = ones(nj-2,ni-2)*gamma;
gamma_w = gamma_e ;
gamma_n = ones(nj-2,ni-2)*gamma;
gamma_s = gamma_n ;

%% Temperature -Matrix
T=zeros(nj,ni);
Tf_x=zeros(nj-2,ni-1); %Temperatures at the faces
Tf_y=zeros(nj-1,ni-2);
for a= 4:nj-1
    if (a<nj-inlet_node)
    T(a,1)=10; %B.C at wall 1
    else
    T(a,1)=T_A; %inlet   
    end
end
for b =1:nj
if (bc_flag == 1)
T(b,ni)= T_B ;
end
end

T(2:nj-1,2:ni-1) = T_i ;   
B=zeros(nj-2,ni-2);
for s= 1:nj-2
    for z= 1:ni-2
 B(s,z)=  s_bar*4*xn(z)*yn(s);
    end
end
%%
aE=zeros(nj-2,ni-2);
aW=zeros(nj-2,ni-2);
SP = zeros(nj-2,ni-2) ;
SU = zeros(nj-2,ni-2) ;

for z= 1:nj-2
    for s=1:ni-2
    F_e(z,s)=rho*u_e(z,s);
    F_w(z,s)=rho*u_w(z,s);
    D_e(z,s)=gamma_e(z,s)*(2*yn(z))/Xpe_delta(s);
    D_w(z,s)=gamma_w(z,s)*(2*yn(z))/Xwp_delta(s);
    end
end

for z= 1:nj-2
     for s=1:ni-2
aW_1= D_w(z,s)+fwx(s)*F_w(z,s);
aW_2 = D_w(z,s)+F_w(z,s) ;
aw(1) = aW_1 ;
aw(2) = aW_2 ;
aw(3) = D_w(z,s) ;
 aW(z,s) = max (aw);
 if (s==1)
 aW(z,s)= 0 ;
 SP(z,s) = - max (aw)+ SP(z,s);
 SU(z,s) = B(z,s) - SP(z,s)*T(z+1,s);
 end
     end
end

for z= 1:nj-2
     for s=1:ni-2
aE_1=  D_e(z,s) -fex(s)*F_e(z,s);
aE_2 = D_e(z,s) - F_e(z,s);
ae(1) = aE_1 ;
ae(2) = aE_2 ;
ae(3) = D_e(z,s) ;
 aE(z,s) = max(ae);
 if (s==ni-2)
 aE(z,s)= 0 ;
 SP(z,s) = - max(ae)+ SP(z,s);
 SU(z,s) = B(z,s) - SP(z,s)*T(z+1,s+2);
 end
     end
end

aW(1:outlet_node,1) = 0;
SP(1:outlet_node,1) = 0;
if (bc_flag == 0)
SU(:,ni-2) = 0 ;
aE(:,ni-2) = 0 ;
SP(:,ni-2) = 0;
end
%%
aN=zeros(nj-2,ni-2);
aS=zeros(nj-2,ni-2);

for z= 1:nj-2
    for s=1:ni-2
    F_n(z,s)=rho*v_n(z,s);
    F_s(z,s)=rho*v_s(z,s);
    D_n(z,s)=gamma_n(z,s).*(2*xn(s))./Ypn_delta(z);
    D_s(z,s)=gamma_s(z,s).*(2*xn(s))./Ysp_delta(z);
    end
end

for z= 1:nj-2
     for s=1:ni-2
aS_1= D_s(z,s)+fsy(s)*F_s(z,s);
aS_2 = D_s(z,s)+F_s(z,s);
as(1) = aS_1 ;
as(2) = aS_2 ;
as(3) = D_s(z,s) ;
aS(z,s) = max (as);
     end
end

for z= 1:nj-2
     for s=1:ni-2
aN_1=  D_n(z,s)-fny(s)*F_n(z,s);
aN_2 = D_n(z,s) - F_n(z,s);
an(1) = aN_1 ;
an(2) = aN_2 ;
an(3) = D_n(z,s) ;
aN(z,s) = max (an);
     end
end

aS(1,:) = 0 ;
aN(nj-2,:) = 0 ;

%%
aP = aE+aW+aN+aS-SP;


%%
grad_x = zeros(nj-2,ni-2);
grad_y = zeros(nj-2,ni-2);
h_A = y(nj-1)-y(nj-4);
iter=0;
[xnc_m,ync_m] = meshgrid(xnc(2:(ni-1)),ync(2:(nj-1)));
cnv=0.1;
if (s_flag==1)
 while (cnv>error)
a_i = aN;
b_i = aS;
D_i = aP;
for s= 1:nj-2

       for z=1:ni-2
C_i(s,z) = aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+SU(s,z);
       end
end

for z = 1:ni-2
    for p = 1:ni-2
        if (p==1)
        A_i(p,z) = a_i(p,z)./(D_i(p,z)-(b_i(p,z)*1));
        C_i_tilde(p,z) = (C_i(p,z)+(1*b_i(p,z)))./(D_i(p,z)-(b_i(p,z)*1));
        else
    A_i(p,z) = a_i(p,z)./(D_i(p,z)-(b_i(p,z)*A_i(p-1,z)));
    C_i_tilde(p,z) = (C_i(p,z)+(C_i_tilde(p-1,z)*b_i(p,z)))./(D_i(p,z)-(b_i(p,z)*A_i(p-1,z)));
        end
    end
end

for i= 1:ni-2
    for j=1:nj-2
  T(nj-j,i+1) = A_i(nj-1-j,i)* T(nj+1-j,i+1) +C_i_tilde(nj-1-j,i) ; 
    end
end

T(nj,:)=T(nj-1,:); %B.C 
T(1,:)=T(2,:); %B.C 
if (bc_flag == 0)
T(:,ni)=T(:,ni-1); %B.C east wall
end
T(1:3,1)=T(1:3,2); %outlet B.C

R=0;
for s= 1:nj-2
        for z=1:ni-2
         R=sum(sum(abs((aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+SU(s,z))-(aP(s,z)*T(s+1,z+1)))))+R;
        end
end

iter=iter+1;
F = rho*u_w(nj-2)*h_A*(T_A-T(2,1));
conv_hist(iter)=R/F;
cnv=conv_hist(iter);

end
elseif(s_flag==0)
 while (cnv>error)
% for i=1:iter
    for s= 1:nj-2
        for z=1:ni-2
            T(s+1,z+1)=(aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+SU(s,z))/(aP(s,z));                        
        end
    end

T(nj,:)=T(nj-1,:); %B.C 
T(1,:)=T(2,:); %B.C 
if (bc_flag == 0)
T(:,ni)=T(:,ni-1); %B.C east wall
end
T(1:3,1)=T(1:3,2); %outlet B.C
% end
R=0;
for s= 1:nj-2
        for z=1:ni-2
         R=sum(sum(abs((aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+SU(s,z))-(aP(s,z)*T(s+1,z+1)))))+R;
        end
end

iter=iter+1;
F = rho*u_w(nj-2)*h_A*(T_A-T(2,1));
conv_hist(iter)=R/F;
cnv=conv_hist(iter);

 end

%global conservation
Tf_x = zeros(nj-2,ni-1) ;
for s= 1:nj-2
        for z=1:ni-3
Tf_x(s,z+1)=interp1(xnc,T(s+1,:),x(z+1)) ;
        end
end
Tf_x(:,1)= T(2:nj-1,1);
Tf_x(:,ni-1)= T(2:ni-1,ni);

Tf_y = zeros(nj-1,ni-2) ;
for z= 1:ni-2
        for s=1:nj-3
 Tf_y(s+1,z)=interp1(ync,T(:,z+1),y(s+1)) ;
        end
end
Tf_y(1,:)= T(1,2:ni-1);
Tf_y(nj-1,:)= T(nj,2:ni-1);

for s= 1:nj-2
for z=1:ni-2
         grad_x(s,z) = -gamma_node(s,z)*(Tf_x(s,z+1)-Tf_x(s,z))./(x(z+1)-x(z));
         grad_y(z,s) = -gamma_node(s,z)*(Tf_y(z+1,s)-Tf_y(z,s))./(y(z+1)-y(z));
end
end

glb_cnv_x = sum( sum(grad_x(:,1))+sum(grad_x(:,ni-2))+sum(grad_x(1,:))+sum(grad_x(nj-2,:))) ;
glb_cnv_y = sum(sum(grad_y(:,1))+sum(grad_y(:,ni-2))+sum(grad_y(1,:))+sum(grad_y(nj-2,:)));

end
figure(1)
hold on
T_node = T(2:(nj-1),2:(ni-1));
contourf(xnc_m,ync_m,T_node)
vec= 10;
quiver(xnc,ync,u2d,v2d,vec)
% xlim([1.5,3])
% ylim([1.2,2])
% figure (2)
% contourc(T)
% pcolor(xnc_m,ync_m,T_node);
% hold on;
% shading interp;
% contour(xnc_m,ync_m,T_node,'LineColor','k')
% contourf(T)
% title('Temperature') ;
%     xlabel({'x direction [m] '});
%     ylabel({'y direction [m]'});
% figure(2)
% plot(grad_x(:,1),ync_m,'r',zeros(nj-2,1),ync_m,'b--')
% title('Heat Flux along the left wall(Dirichlet B.C)') ;
%     xlabel({' Heat Flux '});
%     ylabel({'y direction [m]'});
% legend('Heat Flux in x direction','location','southwest')    
% ylim([0,0.2])
% xlim([-3.5,0.5])
% figure(3)
% plot(grad_y(:,1),ync_m,'r',zeros(nj-2,1),ync_m,'b--')
% % quiver(xnc_m,ync_m,grad_x,grad_y)
% title('Heat Flux along the left wall(Dirichlet B.C)') ;
%     xlabel({' Heat Flux '});
%     ylabel({'y direction [m]'});
% legend('Heat Flux in y direction','location','southwest')
% % ylim([0,2])
toc