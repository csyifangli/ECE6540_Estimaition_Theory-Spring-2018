% SVD is not implemented in this code
% SVD implementation, please see task_SVD 
%========================MODEL SETUP=========================
clc;
clear all;
format long 
doa=[15 20]/180*pi; %actual DOA
N=10;%# Snapshots
w=[pi/3 pi/5]';%Frequency of the envelope signal, slow changing
P=length(w); %# of sources
M=10; %# antennas
lambda=150;%Wavelength
d=lambda/2;%half-wave spacing
snr=20;%SNR
D=zeros(M,P); 
gama=60 ;% penalty factor: tune this parameter to obtain the best DOA result


%-------------------construct the snapshots-------------------
for k=1:P
D(:,k)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]'); %Assignment matrix
end
 xx=2*exp(j*(w*[1:N]));% N=1-> single snapshot
x_s=D*xx;
x_s=x_s+awgn(x_s,snr);%nosie-corrupted sensed signals
%===========================MUSIC algorithm============================
R=x_s*x_s'; %Data covarivance matrix
[U,V]=eig(R); %Find the eigenvalues and eigenvectors of R
NN=U(:,1:M-P); %Estimate noise subspace
theta=-90:0.5:90; %Peak search
for ii=1:length(theta)
SS=zeros(1,length(M));
for jj=0:M-1
SS(1+jj)=exp(j*2*jj*pi*d*sin(theta(ii)/180*pi)/lambda);
end
PP=SS*NN*NN'*SS';
Pmusic(ii)=abs(1/PP);
end
Pmusic=10*log10(Pmusic/max(Pmusic)); %Spatial spectrum function


% =================== L-1 SVD algorithm=============================
%----------------optimization parameters-----------------
Theta_grid=-90:90;
N_grid=size(Theta_grid,2);%  # of grid used in sparse-OPT
D_ext=zeros(M,N_grid);% grid-based steering matrix

Theta_grid=Theta_grid/180*pi;

for k=1:N_grid
D_ext(:,k)=exp(-j*2*pi*d*sin(Theta_grid(k))/lambda*[0:M-1]'); 
end

A = -D_ext;
A_hat=A;
for i=1:N-1
 A_hat=blkdiag(A_hat,A); % stich the A s to get A_hat 
end


% this three-dimensional matrix causes memory overflow
% find a way to express it as 2-dimensional matrices
% this also validified the paper's assertion that SVD is necessary 
%to reduce the computational complecity

M_hat=zeros(N*N_grid,N*N_grid,N_grid);
for i=1:N_grid
    for j=1:N
        M_hat(i+(j-1)*N_grid,i+(j-1)*N_grid,i)=1;
    end
end

y=zeros(M*N,1);
for i=1:M
    for j=N
        y(i+(j-1)*M)=x_s(i,j);
    end
end
Id=ones(N_grid,1);
%-------------------------CVX SOCP--------------------------------
    cvx_begin
        variable s(N*N_grid);
        variable r(N_grid);
        variable p;
        variable q;
        minimize(p+gama*q);
        subject to
            norm(A_hat*s+y)<=sqrt(p); %||A_hat*s+y||2<=sqrt(p)
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
for j=1:N_grid
    P_sparse_dB(j)=10*log10( P_sparse(j)/max( P_sparse));
end
Grid_index=-90:90;
%===================plot results======================================
plot(theta,Pmusic,'-k',Grid_index,P_sparse_dB,'r');
%plot(theta,Pmusic,'-k' );
xlabel('angle \theta/degree')
ylabel('spectrum function P(\theta) /dB')
title('DOA estimation l-1 SVD(red) VS MUSIC(black)')
grid on




