function [ber, numBits] = exp1_MonteCarlo_With_Interleaving(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;

% Initialize variables related to exit criteria.
berVec = zeros(3,1); % Updated BER values


n = 7;
k = 4;
M = 2;
K = log2(M);
Len = 1e4;
alphabet = 0 : M - 1; 
A = 1; 
teta0 = 0; 
Ts = 1/1e5;
Fd = 500;

%g = ones(1,fs)/sqrt(fs);
% --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((berVec(2) < maxNumErrs) && (berVec(3) < maxNumBits))
        
    
        data = randsrc(1, Len, alphabet); 
        enc = encode(data, n,k,'hamming/binary');
%         S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 
%         ss = S_PSK(enc + 1); 
        ss = 2*enc - 1;
          
        int_1=reshape(ss,[],size(ss,1));
        int_2=transpose(int_1);
        int_3=reshape(int_2,1,[]);
        
        sym = int_3*sqrt(k/n);
        
        y = awgn(sym, EbNo + 3);
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
    

      
        berVec(2) = berVec(2) + number;
        berVec(3) = berVec(3) + Len;
   % Check if the user clicked the Stop button of BERTool.
   if (BERTool.getSimulationStop)
      break;
   end

   % --- Proceed with simulation.
   % --- Be sure to update totErr and numBits.
   % --- INSERT YOUR CODE HERE.
end % End of loop

berVec(1) = berVec(2)/berVec(3);

% Assign values to the output variables.
ber = ratio;
numBits = berVec(3);
