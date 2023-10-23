%A
clc
close all
clear

EbN0 = 10;
n = 0;

trainlen = 1e3;
datalen = 1e4;

b_train = randsrc(1,trainlen,[0,1]); 
t_train = 2 * b_train - 1; 

b_data = randsrc(1,datalen,[0,1]); 
t_data = 2 * b_data - 1; 

pre_seq = [b_train b_data];
sequence = [t_train t_data];

channel_A = [0.04 0.05 0.07 0.21 0.5 0.72 0.36 0 0.21 0.03 0.07];

recieved_seq = conv(sequence, channel_A);

r = awgn(recieved_seq, EbN0 + 3);

r_sign = r>0;

real_recieved = r_sign(6 : length(r_sign) - 5);

n = nnz(real_recieved - pre_seq);

BER = n/length(pre_seq);

figure;
hist(r);
%%
%E
spectra_A = fft(channel_A, 1e5);

real_spectra_A = fftshift(abs(spectra_A).^2);

figure;
plot(real_spectra_A);
title("PSD");

figure;
freqz(channel_A);
%%
%B
clc
close all

L = 11;
N = (L-1)/2;
delta = 0.01;
r_real = r(6 : length(r) - 5);
Yk = zeros(1,length(r_real)-10);
coefficient = zeros(length(r_real)-10,L );
coefficient(1,:) = eps*ones(1,L);

for i = 1 : trainlen+datalen -10
    
    Yk(i) = sum(coefficient(i,:).*r_real(i : i+10));

    if i <= trainlen
       
        err = pre_seq(i) - Yk(i);
        zz = delta*err*r_real(i : i+10);
        coefficient(i+1 , :) =  coefficient(i,:)+ zz ;
        
    end
    
    err = sign(Yk(i)) - Yk(i);
    zz = delta*err*r_real(i : i+10);
    coefficient(i+1 , :) =  coefficient(i,:)+ zz; 
   
end


% trainlen = 1e2;
% datalen = 1e4;
% 
% b_train = randsrc(1,trainlen,[0,1]); 
% t_train = 2 * b_train - 1; 
% 
% b_data = randsrc(1,datalen,[0,1]); 
% t_data = 2 * b_data - 1; 
% 
% pre_seq = [b_train b_data];
% sequence = [t_train t_data];
% 
% channel_A = [0.04 0.05 0.07 0.21 0.5 0.72 0.36 0 0.21 0.03 0.07];


equalized_signal = conv(r, coefficient(length(coefficient),:));

r_sign = equalized_signal>0;

real_recieved = r_sign(11 : length(r_sign) - 10);

n = nnz(real_recieved - pre_seq);

BER = n/length(pre_seq);

figure;
hist(r);
