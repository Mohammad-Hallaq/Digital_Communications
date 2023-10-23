function S = creatSquareConstellation(M)
% Written by: Majdi Msallam
if M == 4
S = [-1+1i -1-1i 1+1i 1-1i];
elseif M == 8
S = [-3+1i -3-1i -1+1i -1-1i 1+1i 1-1i 3+1i 3-1i];
else
K = log2(M); % number of bits per symbol
if mod(K, 2) == 0 % perfect square constellation M = 64, 256, 1024, ...
m = sqrt(M); % number of projections on each axis
x = -m + 1 : 2 : m - 1; % real part
S = [];
for ii = 1 : m
% start at the left-down corner of constellation
temp = [S (x(ii) + 1i * x)];
S = temp;
end
else% square constellation with gaps at corners: M = 32, 128, 512, ...
s1 = log2(M) - 1; % inner square dimension
e = (M - 2 ^ s1) / 4 / sqrt(2 ^ s1); % extra points dimension
m = sqrt(2 ^ s1) + 2 * e; % outer square dimension
x = -m + 1 : 2 : m - 1; % projections on the real axis
S = [];
% start at the left-down corner of constellation
for ii = 1 : e % left part of constellation
temp = [S (x(ii) + 1i * x(e + 1 : m - e))]; S = temp;
end
for ii = e + 1 : m - e % middle part
temp = [S (x(ii) + 1i * x)]; S = temp;
end
for ii = m - e + 1 : m % right part
temp = [S (x(ii) + 1i * x(e + 1 : m - e))]; S = temp;
end
end
end
p_qam = sum(abs(S) .^ 2) / M;
S = S / sqrt(p_qam); % scaling to unity mean power