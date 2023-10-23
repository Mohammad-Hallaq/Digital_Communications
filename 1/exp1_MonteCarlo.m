function [ber, numBits] = exp1_MonteCarlo(EbNo, maxNumErrs, maxNumBits)
% Import Java class for BERTool.
import com.mathworks.toolbox.comm.BERTool;

% Initialize variables related to exit criteria.
berVec = zeros(3,1); % Updated BER values
Len = 1e3;
alphabet = 0 : M - 1; 

A = 1; 
teta0 = 0; 
S_PSK = A * exp(1i * (2 * pi * alphabet / M + teta0)); 

% --- Set up parameters. ---
% --- INSERT YOUR CODE HERE.
% Simulate until number of errors exceeds maxNumErrs
% or number of bits processed exceeds maxNumBits.
while((berVec(2) < maxNumErrs) && (berVec(3) < maxNumBits))
    
        data = randsrc(1, Len, alphabet);% information bits
        s = S_PSK(data + 1); 
        r = awgn(s, EbNo + 10*log2(M)); 
        
        for i=1:length(r)
             for j=1:M
        
                  dist_Mat(i,j) = abs(r(i) - S_PSK(j));
        
             end
    
        end
        [dist_min_Mat, index] = min(dist_Mat,[],2);
        recieved_data = index -1;
        
        [number,ratio] = biterr(recieved_data', data, K);
        
         berVec(2) = berVec(2) + number;
         berVec(3) = berVec(3) + Len*log2(M);
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
ber = berVec(1);
numBits = berVec(3);