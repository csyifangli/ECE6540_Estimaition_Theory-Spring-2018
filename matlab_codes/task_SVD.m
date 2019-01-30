% SVD is not implemented in this code
% SVD implementation, please see task_SVD 
%========================MODEL SETUP=========================
clc;
clear all;
format long 
doa=[75 ]/180*pi; %actual DOA
N=100;%# Snapshots
w=[pi/3 ]';%Frequency of the envelope signal, slow changing
P=length(w); %# of sources
M=10; %# antennas
lambda=150;%Wavelength
d=lambda/2;%half-wave spacing
snr=20;%SNR
D=zeros(M,P); 



%-------------------construct the snapshots-------------------
for k=1:P
D(:,k)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]'); %Assignment matrix
end
 xx=2*exp(j*(w*[1:N]));% N=1-> single snapshot
x_s=D*xx;
x_s=x_s+awgn(x_s,snr);%nosie-corrupted sensed signal matrix including multiple snapshots

% %===========================MUSIC algorithm============================
R=x_s*x_s'; %Data covarivance matrix
[G,H]=eig(R); %Find the eigenvalues and eigenvectors of R
NN=G(:,1:M-P); %Estimate noise subspace && MUSIC needs to know the # of sources P
theta=-90:90; %Peak search
for ii=1:length(theta)
SS=zeros(1,length(M));
for jj=0:M-1
SS(1+jj)=exp(j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP=SS*NN*NN'*SS';
Pmusic(ii)=abs(1/PP);
end
Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function
plot(theta,Pmusic,'-k');
hold on;

% =================== L-1 SVD algorithm=============================
%--------------------SVD transformation-------------------------------
[U,L,V]=svd(x_s); %  SVD decomposition, L is diagonal
K=rank(x_s,0.1*L(1,1));% return the # of singular values that are greater than 0.1*L(1,1)
% signal space has the largest K singular values and noise space has very
% small values. K should be equal to P(# of sources)
D_K=zeros(N,K);
for i=1:K
    D_K(i,i)=1;
end
Y_svd=U*L*D_K;% noise matrix trans is unnecessary cuz unitary tansform doesm't change the distri. of Gaussain noise
D_K
%----------------optimization parameters-----------------
N_grid=181;%  # of grid used in sparse-OPT
D_ext=zeros(M,N_grid);% grid-based steering matrix
Theta_grid=-90:90;
Theta_grid=Theta_grid/180*pi;

for k=1:N_grid
D_ext(:,k)=exp(-j*2*pi*d*sin(Theta_grid(k))/lambda*[0:M-1]'); 
end

A = -D_ext;
A_hat=A;
for i=1:K-1
 A_hat=blkdiag(A_hat,A); % stich the A s to get A_hat 
end

M_hat=zeros(K*N_grid,K*N_grid,N_grid);
for i=1:N_grid
    for j=1:K
        M_hat(i+(j-1)*N_grid,i+(j-1)*N_grid,i)=1;
    end
end
y=zeros(M*K,1);
for i=1:M
    for j=K
        y(i+(j-1)*M)=Y_svd(i,j);
    end
end
Id=ones(N_grid,1);
for gama=700 % find the best penalty factor
%-------------------------CVX SOCP--------------------------------
    cvx_begin
        variable s(K*N_grid);
        variable r(N_grid);
        variable p;
        variable q;
        minimize(p+gama*q);
        subject to
            norm(A_hat*s+y)<=sqrt(p); %||A_hat*s+y||2<=sqart(p)
            Id'*r<=q;  %sum(r)<=q
            for i=1:N_grid
                norm(M_hat(:,:,i)*s)<=r(i);%||M_hat*s||2<=r_i
            end        
    cvx_end
    
%-----------------------------------------------------------
P_sparse=[];% Sparse- based spectrum  
P_sparse_dB=[];
for j=1:N_grid
    P_sparse(j)=norm(M_hat(:,:,j)*s);
end
P_sparse_dB=10*log10( P_sparse/max( P_sparse));
Grid_index=-90:90;
plot(Grid_index,P_sparse_dB,'r');
hold on
end
%===================plot results======================================
% plot(theta,Pmusic,'-k',Grid_index,P_sparse_dB,'r');
 xlabel('Angle \theta /degree')
ylabel('Spatial Spectrum P(\theta) /dB')
% title('DOA estimation l-1 SVD(red) VS MUSIC(black)')
% grid on




