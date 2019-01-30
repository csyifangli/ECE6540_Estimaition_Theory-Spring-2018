


% ==================single snapshot L-1 SVD DOA algorithm======================
clc;
clear all;
format long %The data show that as long shaping scientific
doa=[-4 4]/180*pi; %Direction of arrival
N=1;%# Snapshots
w=[pi/3 pi/4]';%Frequency of the envelope signal, slow changing
P=length(w); %The number of signal
M=10;%Number of array elements
lambda=150;%Wavelength
d=lambda/2;%half-wave spacing
snr=20;%SNR
D=zeros(M,P); %To creat a matrix with P row and M column

for k=1:P
D(:,k)=exp(-j*2*pi*d*sin(doa(k))/lambda*[0:M-1]'); %Assignment matrix
end

D