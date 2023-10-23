clc
close all
clear

M = 4;
K = log2(M);

Len = 1e3;

alphabet = 0 : M - 1; 
data = randsrc(1, Len, alphabet); 

SNR = 7;

% PSK
A = 1; 
teta0 = 0; 
S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
s = S_PSK(data + 1); 
r = awgn(s, SNR); 
figure(1); scatter(real(r), imag(r)); grid on; hold on
title(['Constellation Points for ' num2str(M) '-PSK Modulation'])
scatter(real(S_PSK), imag(S_PSK), '*r');
legend('Received Symbols', 'Constellation Points');

for i=1:M
    for j=1:M
        if i~=j
          d(i,j) = sqrt((real(S_PSK(i)) - real(S_PSK(j)))^2 + (imag(S_PSK(i)) - imag(S_PSK(j)))^2);
        end
        if i == j
            d(i,j)=1e5;
        end
    end
end

dmin = min(min(d));



for i=1:length(r)
    for j=1:M
        
        dist_Mat(i,j) = abs(r(i) - S_PSK(j));
        
    end
    
end

[dist_min_Mat, index] = min(dist_Mat,[],2);

recieved_data = index -1;

[number,ratio] = biterr(recieved_data', data, K);
%%
% QAM
clc

M = 4;
K = log2(M);

alphabet = 0 : M - 1; 
data = randsrc(1, Len, alphabet); 

SNR = 30;

S_QAM = creatSquareConstellation(M);
s = S_QAM(data + 1);
r = awgn(s, SNR);
figure(2); scatter(real(r), imag(r)); grid on; hold on
title(['Square Constellation for ' num2str(M) '-QAM Modulation'])
scatter(real(S_QAM), imag(S_QAM), '*r');
legend('Received Symbols', 'Constellation Points');

for i=1:M
    for j=1:M
        if i~=j
          d(i,j) = sqrt((real(S_QAM(i)) - real(S_QAM(j)))^2 + (imag(S_QAM(i)) - imag(S_QAM(j)))^2);
        end
        if i == j
            d(i,j)=1e5;
        end
    end
end

dmin = min(min(d));

for i=1:length(r)
    for j=1:M
        
        dist_Mat(i,j) = abs(r(i) - S_QAM(j));
        
    end
    
end

[dist_min_Mat, index] = min(dist_Mat,[],2);

recieved_data = index -1;


%%

