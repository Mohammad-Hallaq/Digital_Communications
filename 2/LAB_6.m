clc
clear 
close all

[H,G,n,k] = hammgen(3);

msg = randsrc(1,k,[0 1]);

code_word = mod(msg*G, 2);

S = mod(code_word*H', 2);

error_code = mod(code_word + [0 0 0 1 0 0 0], 2);

S2 = mod(error_code*H', 2);

T = syndtable(H);

%%
clc

M = 2;
K = log2(M);

Ts = 1/1e5;
Fd = 500;
Len = 1e3;

alphabet = 0 : M - 1; 

data = randsrc(1, Len, alphabet); 

enc = encode(data, 7,4,'hamming/binary');

%interleaver
dd = reshape(enc,7,[]);

pad = [dd zeros(7,2)];

for i=1:length(pad(:,1))
    
    ddd(:,(1:36)+(i-1)*(length(pad(1,:)))/7) = reshape(pad(i,:),7,[]);
    
end

SNR = 3;

% PSK
A = 1; 
teta0 = 0; 
S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
ss = S_PSK(ddd + 1); 

sym = ss*sqrt(k/n);
series = reshape(sym,1764,[]);
ch = rayleighchan(Ts, Fd);
ch2 = abs(filter(ch,ones(1,length(series))));
y = series.*ch2;
r = awgn(series', 0 + 3);
        


        

for i=1:length(r)
    for j=1:M
        
        dist_Mat(i,j) = abs(r(i) - S_PSK(j));
        
    end
    
end

[dist_min_Mat, index] = min(dist_Mat,[],2);
recieved_data = index -1;


de_one = reshape(recieved_data,7,[]);

for count5 = 1:length(de_one(:,1))

for count6=1:ceil(250/7)
 
    rev_pad_samp(count5,((count6-1)*7)+(1:7)) = de_one (1:7, count6+ (count5 - 1)*36 );
    
end

end

rev_samp = rev_pad_samp(:,1:length(rev_pad_samp)-2);

 for num =1:length(rev_samp)
     
     de_inter((num-1)*7 + (1:7)) = rev_samp(1:7,num);
     
 end

 
ch = de_inter-enc;
 
dec = decode(de_inter, n,k, 'hamming/binary');

[number,ratio] = biterr(dec, data, K);

SNR = [-8:0.5:8];
for jj =1:length(SNR)
    
    S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
    ss = S_PSK(enc + 1); 

    sym = ss*sqrt(k/n);

    r = awgn(sym, SNR(jj)); 


for i=1:length(r)
    for j=1:M
        
        dist_Mat(i,j) = abs(r(i) - S_PSK(j));
        
    end
    
end

[dist_min_Mat, index] = min(dist_Mat,[],2);
recieved_data = index -1;

dec = decode(recieved_data, n,k, 'hamming/binary');

[number(jj),ratio(jj)] = biterr(dec', data, K);
end

plot(ratio);
grid on ;
%%
clear
clc

data = randsrc(1, 1000, [0,1]); 

Ts = 1/1e5;
Fd = 500;

ch = rayleighchan(Ts, Fd,[0 1e-5 5*1e-5], [1 0.7 0.1]);

ch.storeHISTORY = 1;

y = filter(ch, data);

figure;
plot(ch);
%%
clc

Ts = 1/1e5;
Fd = 500;

ch = rayleighchan(Ts, Fd);

ch.storeHISTORY = 1;

y = filter(ch, sym);

SNR =9;

r = awgn(y, SNR); 

for i=1:length(r)
    for j=1:M
        
        dist_Mat(i,j) = abs(r(i) - S_PSK(j));
        
    end
    
end

[dist_min_Mat, index] = min(dist_Mat,[],2);
recieved_data = index -1;

dec = decode(recieved_data, n,k, 'hamming/binary');

[number,ratio] = biterr(dec', data, K);
%%
clear all
clc
%data = randsrc(1e3,1,[0,1]);

count2 = 1;
for count1=1:length(enc)
    
    if mod(count1,7) == 0
        
        samp(7,count2) = enc(count1);
       
        count2 = count2 + 1;
               
    else
        
        samp(mod(count1,7),count2) = enc(count1);
       
    end
    
end

pad_samp = [samp zeros(7,2)];

for count3 = 1:length(pad_samp(:,1))

for count4=1:ceil(250/7)
 
    inter (1:7, count4+ (count3 - 1)*36 ) = pad_samp(count3,((count4-1)*7)+(1:7)); 
    
end

end


%%
% PSK
A = 1; 
teta0 = 0; 
alphabet = 0 : M - 1; 
S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
ss = S_PSK(inter + 1); 
sym = ss*sqrt(k/n);

for ll=1:7
for ii=1:length(sym)
    for kj=1:M
        
        dist_Mat(kj,ll, ii) = abs(sym(ll,ii) - S_PSK(kj));
        
    end
    
end
end

[dist_min_Mat, index] = min(dist_Mat,[],1);
recieved_data = index -1;

rr = reshape(recieved_data, [7, 252]);

for count5 = 1:length(rr(:,1))

for count6=1:ceil(250/7)
 
    rev_pad_samp(count5,((count6-1)*7)+(1:7)) = rr (1:7, count6+ (count5 - 1)*36 );
    
end

end

rev_samp = rev_pad_samp(:,1:length(rev_pad_samp)-2);

 for num =1:length(rev_samp)
     
     dennc((num-1)*7 + (1:7)) = rev_samp(1:7,num);
     
 end
 %%
 
 
n = 7;
k = 4;
M = 2;
K = log2(M);
Len = 1e3;
alphabet = 0 : M - 1; 
A = 1; 
teta0 = 0; 
Ts = 1/1e5;
Fd = 500;
data = randsrc(1, Len, alphabet); 
        enc = encode(data, n,k,'hamming/binary');
        dd = reshape(enc,7,[]);
        pad = [dd zeros(7,2)];

        for i=1:length(pad(:,1))

            inter(:,(1:36)+(i-1)*(length(pad(1,:)))/7) = reshape(pad(i,:),7,[]);

        end
        
        S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
        ss = S_PSK(inter + 1); 
        sym = ss*sqrt(k/n);
        series = reshape(sym,1764,[]);
        
        y = awgn(series, EbNo + 3);
        ch = rayleighchan(Ts, Fd);
        ch2 = abs(filter(ch,ones(1,length(y))));
        r = y.*ch2;
        
        
        for i=1:length(r)
            for j=1:M

                dist_Mat(i,j) = abs(r(i) - S_PSK(j));

            end

        end

        [dist_min_Mat, index] = min(dist_Mat,[],2);
        recieved_data = index -1;
        
        de_one = reshape(recieved_data,7,[]);

        for count5 = 1:length(de_one(:,1))

        for count6=1:ceil(250/7)

            rev_pad_samp(count5,((count6-1)*7)+(1:7)) = de_one (1:7, count6+ (count5 - 1)*36 );

        end

        end

        rev_samp = rev_pad_samp(:,1:length(rev_pad_samp)-2);

         for num =1:length(rev_samp)

             de_inter((num-1)*7 + (1:7)) = rev_samp(1:7,num);

         end
           
        
        dec = decode(de_inter, n,k, 'hamming/binary');
        [number,ratio] = biterr(dec, data, K);
      
        
%%

n = 7;
k = 4;
M = 2;
K = log2(M);
Len = 1e3;
alphabet = 0 : M - 1; 
A = 1; 
teta0 = 0; 
Ts = 1/1e5;
Fd = 500;
  data = randsrc(1, Len, alphabet); 
        enc = encode(data, n,k,'hamming/binary');
%         S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
%         ss = S_PSK(enc + 1); 
        ss = 2*enc - 1;
          
        int_1=reshape(ss,[],size(ss,1));
        int_2=transpose(int_1);
        int_3=reshape(int_2,1,[]);
        
        sym = int_3*sqrt(k/n);
        
        y = awgn(sym, 0 + 3);
        ch = rayleighchan(Ts, Fd);
        ch2 = abs(filter(ch,ones(1,length(y))));
        r = y.*ch2;
        
        
%         for i=1:length(r)
%             for j=1:M
% 
%                 dist_Mat(i,j) = abs(r(i) - S_PSK(j));
% 
%             end
% 
%         end

%         [dist_min_Mat, index] = min(dist_Mat,[],2);
%         recieved_data = index -1;
        recieved_data = (sign(r)+1)/2;
        
        de_int_1=reshape(recieved_data,[],size(recieved_data,1));
        de_int_2=transpose(de_int_1);
        de_int_3=reshape(de_int_2,1,[]);
        
        dec = decode(de_int_3, n,k, 'hamming/binary');
        [number,ratio] = biterr(dec, data, K);
        
        
        
        %%
        m = randsrc(1,1000,[0 1]);
       es = encode_func(randsrc(1,1000,[0 1]));
       
       ees = encode_func(es);
    