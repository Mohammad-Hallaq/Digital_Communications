function [ series ] = decode_func( input_vec)

[H,G,n,k] = hammgen(3);

trt = syndtable(H); 
t=length(input_vec)/n;

series=0;
for i=1:t
 code_word=input_vec(1+(i-1)*n:i*n);
 syndrome = rem( code_word* H',2);
 syndrome_de = bi2de(syndrome,'left-msb'); 
 corrvect = trt(1+syndrome_de,:); 
 correctedcode = rem(corrvect+code_word,2);
 series=[series correctedcode(n-k+1:n)];
end
series=series(2:length(series));
end

