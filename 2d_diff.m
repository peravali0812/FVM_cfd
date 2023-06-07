%------------------------------------------------------------------------------------------
%                      2D DIFFSION EQUATION USING MATLAB(Surya Kiran Peravali)
%                      CHALMERS UNIVERSITY OF TECHNOLOGY, Gothenburg,2014
%------------------------------------------------------------------------------------------

clc
clear all
close all
tic
%Inputs
L =1.5; %length od domian
H = 0.5; %Height of the domain
s_bar = 0; % source term
kfactor=1; %for increasing or decreasing k value
iter = 10000 ; % number of iterations
kx_change=[0.7 1.1];
ky_change=[0.3 0.4];
%----- boundary condition type flag ----%
%-- 0 for Neumann B.C --&
%-- 1 for Dirichlet B.C  --&
bc_flag = 0 ;
mesh_flag = 3 ; %type flag number here1
%----- mesh type flag ----%
%-- 1 for 10x10 equidistant mesh --&
%-- 2 for 20x20 equidistant mesh --&
%-- 3 for 40x40 equidistant mesh --&
%-- 4 for 10x10 refined mesh --&
%-- 5 for 20x20 refined mesh --&    
%-- 6 for 40x40 refined mesh variant 1 --&
%-- 7 for 40x40 refined mesh variant 2 --&
%-- 8 for manually defined size, equidistant mesh --&
ni=100; %no of nodes in x direction
nj=100; %no of nodes in y direction
%% selected mesh type parameters
if mesh_flag == 1
    ni=10; %no of nodes in x direction
    nj=10; %no of nodes in y direction
elseif mesh_flag == 2
    ni=20; %no of nodes in x direction
    nj=20; %no of nodes in y direction
elseif mesh_flag == 3
    ni=40; %no of nodes in x direction
    nj=40; %no of nodes in y direction
elseif mesh_flag == 4    
    load mesh1.mat %10x10 mesh
    x = mesh1(:,1)'./1.5*L; %nondimensionalized
    y = mesh1(:,2)'./0.5*H; %nondimensionalized
elseif mesh_flag == 5    
    load mesh2.mat %20x20 mesh
    x = mesh2(:,1)'/1.5*L; %nondimensionalized
    y = mesh2(:,2)'/0.5*H; %nondimensionalized
elseif mesh_flag == 6    
    load mesh3.mat %40x40 mesh
    x = mesh3(:,1)'/1.5*L; %nondimensionalized
    y = mesh3(:,2)'/0.5*H; %nondimensionalized
elseif mesh_flag == 7    
    load mesh4.mat %40x40 mesh
    x = mesh4(:,1)'/1.5*L; %nondimensionalized
    y = mesh4(:,2)'/0.5*H; %nondimensionalized
end

if ((mesh_flag == 1) ||(mesh_flag == 2)||(mesh_flag == 3)||(mesh_flag == 8))
    n=ni+1;
    m=nj+1;
    x=linspace(0,L,n); % x coordinates for faces
    y = linspace(0,H,m); %y coordinates for faces
elseif ((mesh_flag == 4) ||(mesh_flag == 5)||(mesh_flag == 6)||(mesh_flag == 7))
    n=size(x,2);
    m=size(y,2);
    ni=n-1;
    nj=m-1;
else
    disp('Wrong input for mesh flag')
end
%%
xn =zeros(1,ni);
yn =zeros(1,nj);
xnc =zeros(1,ni);
ync =zeros(1,nj);
for i = 1:(length(x)-1)
    xn(i) = (x(i+1)- x(i))/2; %delta x/2
    xnc(i)= x(i+1)-xn(i); %node coordinate x direction
end
for j = 1:(length(y)-1)
    yn(j) = (y(j+1)-y(j))/2; %delta y/2
    ync(j)= y(j+1) -yn(j); %node coordinate in y direction
end

kx_ratio= [kx_change(1)/L kx_change(2)/L];
ky_ratio= [ky_change(1)/H ky_change(2)/H];

%k at nodes
kn = ones(1,ni)*kfactor; 
knode = zeros(nj,ni);
%%Surya%Burak Case
for z = 1:(length(y)-1)
    for s = 1:(length(x)-1)
if ((x(s)>kx_change(1)) && (x(s)<kx_change(2)) && (y(z)>0.3) &&(y(z)<0.4))
    knode(z,s)=0.01.*kn(s);
else
    knode(z,s)=20.*kn(s);
end
    end
end
%% Behrads case
% knode = zeros(nj,ni);
%  for a=1:ni
%     kNodex(a)=5.*(1+100.*(xnc(a)./L)).*kn(a);
%     knode(:,a)=kNodex(a);
%  end
%%

%k at faces
Xpe_delta = zeros(1,ni);
Xwp_delta = zeros(1,ni);
fex = zeros(1,ni);
fwx = zeros(1,ni);
k_e = zeros(nj,ni);
k_w = zeros(nj,ni); 
for j= 1:(length(ync))
for i = 1:(length(xnc)-1)
    Xpe_delta(i)= xnc(i+1)-xnc(i);
    Xpe_delta(ni)= xn(ni);
    Xwp_delta(1) = xn(1);
    Xwp_delta(length(xnc)+1-i)= xnc(length(xnc)+1-i)-xnc(length(xnc)-i);
    fex= xn./Xpe_delta;
    fwx= xn./Xwp_delta;
    k_e(j,i) = (fex(i)*knode(j,i+1)) + ((1-fex(i))*knode(j,i));
    k_e(j,ni) = (fex(ni)*knode(j,ni)) + ((1-fex(ni))*knode(j,ni));
    k_w(j,1)= (fwx(1)*knode(j,1)) + ((1-fwx(1))*knode(j,1));
    k_w(j,length(xnc)+1-i) = (fwx(length(xnc)+1-i)*knode(j,length(xnc)-i))+...
        ((1-fwx(length(xnc)+1-i))*knode(j,length(xnc)+1-i));
end
end
Ypn_delta = zeros(1,nj);
Ysp_delta = zeros(1,nj);
fny = zeros(1,nj);
fsy = zeros(1,nj);
k_n = zeros(nj,ni);
k_s = zeros(nj,ni);

for i = 1:(length(xnc))
for j = 1:(length(ync)-1)
    Ypn_delta(j)=ync(j+1)-ync(j);
    Ypn_delta(nj)= yn(nj);
    Ysp_delta(1) = yn(1);
    Ysp_delta(length(ync)+1-j)= ync(length(ync)+1-j)-ync(length(ync)-j);
    fny = yn./Ypn_delta;
    fsy = yn./Ysp_delta;
    k_n(j,i) = fny(j)*knode(j+1,i) + (1-fny(j))*knode(j,i);
    k_n(nj,i)= fny(nj)*knode(nj,i) + (1-fny(nj))*knode(nj,i);
    k_s(1,j)= (fsy(1)*knode(1,i)) + ((1-fsy(1))*knode(1,i));
    k_s(length(ync)+1-j,i) = (fsy(length(ync)+1-j)*knode(length(ync)-j,i))+...
        ((1-fsy(length(ync)+1-j))*knode(length(ync)+1-j,i));
end
end
%Coefficients
aE=zeros(nj,ni);
aW=zeros(nj,ni);

for z= 1:nj
    for s=1:ni
aE(z,s) = k_e(z,s).*(2*yn(z))./Xpe_delta(s);
aW(z,s) = k_w(z,s).*(2*yn(z))./Xwp_delta(s);
    end
end

%Coefficients
aN=zeros(nj,ni);
aS=zeros(nj,ni);
for s= 1:ni
    for z=1:nj
aN(z,s) = k_n(z,s).*(2*xn(s))./Ypn_delta(z);
if (bc_flag==0)
if z==nj
    aN(z,:)=0;
end
end
aS(z,s) = k_s(z,s).*(2*xn(s))./Ysp_delta(z);
    end
end
aP = aE+aW+aN+aS;

%T-Matrix
T=zeros(nj+2,ni+2);
Tf_x=zeros(nj,ni+1); %Temperatures at the faces
Tf_y=zeros(nj+1,ni);
for a= 1:ni+2
    T(1,a)=10; %B.C at wall 1
end
Tf_y(1,:)=T(1,2:ni+1);
for b=1:nj+2
    T(b,1)=10; %B.C at wall 4
end
Tf_x(:,1)=T(2:nj+1,1); %x-dir. temp. at faces at wall 4
if (bc_flag==1)
 for a= 1:ni+2
    T(nj+2,a)=10; %B.C at wall 3
 end
 Tf_y(nj+1,:)=T(nj+1,2:ni+1);
end
for b=1:nj+1
    T(b+1,ni+2)= 10+(20*sin(pi*y(b))./H); %B.C at wall 2
end
 Tf_x(:,ni+1)= T(2:nj+1,ni+2); %x-dir. temp. at faces at wall 4
%source terms
B=zeros(nj,ni);
for s= 1:nj
    for z= 1:ni
 B(s,z)=  s_bar*4*xn(z)*yn(s);
    end
end
%%
% cnv=1;
%Iteration
R=0;
for i=1:iter
    Ti=T;
    for s= 1:nj
        for z=1:ni
            T(s+1,z+1)=(aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+B(s,z))/(aP(s,z));                        
        end
    end
end
if (bc_flag==0)
T(nj+2,:)=T(nj+1,:); %B.C at wall 3
    Tf_y(nj+1,:)=T(nj+2,2:ni+1);
end 
    
figure(1)
contourf(T)
for s= 1:nj
        for z=1:ni
         R=sum(sum(abs((aE(s,z)*T(s+1,z+2)+aW(s,z)*T(s+1,z)+...
                aN(s,z)*T(s+2,z+1)+aS(s,z)*T(s,z+1)+B(s,z))-(aP(s,z)*Ti(s+1,z+1)))))+R;
        end
end

%Temperatures at faces between the nodes
for s= 1:nj
        for z=1:ni-1
% Tf_x(s,z+1)=interp1(xnc,T(s+1,2:ni+1),x(z+1)) ;
  Tf_x(s,z+1)= (T(s+1,z+2)-T(s+1,z+1))*(x(z+1)-xnc(z))/(xnc(z+1)-xnc(z))+T(s+1,z+1) ;
        end
end
for z=1:ni
        for s= 1:nj-1
% Tf_y(s+1,z)=interp1(ync,T(2:nj+1,z+1),y(s+1)) ;
  Tf_y(s+1,z)= (T(s+2,z+1)-T(s+1,z+1))*(y(s+1)-ync(s))/(ync(s+1)-ync(s))+T(s+1,z+1) ;
        end
end
grad_x = zeros(nj,ni);
grad_y = zeros(nj,ni);
for s= 1:nj
for z=1:ni
        grad_x(s,z) = -knode(s,z)*(Tf_x(s,z+1)-Tf_x(s,z))./(x(z+1)-x(z));
         grad_y(z,s) = -knode(s,z)*(Tf_y(z+1,s)-Tf_y(z,s))./(y(z+1)-y(z));
end
end
F = sum(sum(abs(grad_x)))+sum(sum(abs(grad_y)));
conv=R/F 
[xnc_m,ync_m] = meshgrid(xnc,ync);
figure(1)
T_node = T(2:ni+1,2:nj+1);
contourf(xnc_m,ync_m,T_node)
title('Temperature') ;
    xlabel({'x direction [m] '});
    ylabel({'y direction [m]'});
figure(2)
quiver(xnc_m,ync_m,grad_x,grad_y)
title('Temperature Flux Vectors') ;
    xlabel({'x direction [m] '});
    ylabel({'y direction [m]'});
    
if ((mesh_flag == 1) ||(mesh_flag == 2)||(mesh_flag == 3)||(mesh_flag == 8))
    figure(3)
    hold on
    plot(x,Tf_x(nj*0.2,:),'b',x,Tf_x(nj*0.4,:),'g',x,Tf_x(nj*0.6,:),'r',x,Tf_x(nj*0.8,:),'k')
    plot([kx_ratio(1)*L kx_ratio(1)*L],[0 max(max(Tf_x))],'-.')
    plot([kx_ratio(2)*L kx_ratio(2)*L],[0 max(max(Tf_x))],'-.')
    title('Temperature change in x direction') ;
    xlabel({'x [m] '});
    ylabel({'Temperature [C^0]'});
    legend('T@0.2*H','T@0.4*H','T@0.6*H','T@0.8*H','k-change','location','northeast')
    figure(4)
    hold on
    plot(Tf_y(:,nj*0.2),y,'b',Tf_y(:,nj*0.4),y,'g',Tf_y(:,nj*0.6),y,'r',Tf_y(:,nj*0.8),y,'k')
    plot([0 max(max(Tf_y))],[ky_ratio(1)*H ky_ratio(1)*H],'-.')
    plot([0 max(max(Tf_y))],[ky_ratio(2)*H ky_ratio(2)*H],'-.')
    title('Temperature change in y direction') ;
    ylabel({'y [m]'});
    xlabel({'Temperature [C^0]'});
    legend('T@0.2*L','T@0.4*L','T@0.6*L','T@0.8*L','k-change','location','northeast')
end
figure(5)
z=ones(nj,ni);
surf(xnc_m,ync_m,z)

toc