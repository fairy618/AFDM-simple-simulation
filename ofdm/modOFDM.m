function [dataOut] = modOFDM(dataIn,numSC,cpLen,ofdmSym)

%--------------------------------------------------------------------------
%
%               Modulates random input data into OFDM symbols
%
%--------------------------------------------------------------------------
% Input Arguments: 
% dataIn                         Input data vector
% numSC                          The number of subcarriers used
% cpLen                          Length of the cyclic prefix
% ofdmSym                        No. of ofdm symbols per subframe
%--------------------------------------------------------------------------
% Function returns: 
% dataOut                        Output OFDM symbols
%--------------------------------------------------------------------------

% Calculate variables
cyclicPrefix_start  = numSC - cpLen;

% Perform IFFT
ifftSubcarrier = ifft(dataIn,[],2); 

%Finding cyclic prefix for each subcarrier
for i=1:cpLen
    for j=1:ofdmSym                   
        cyclicPrefix_data(i,j) = ifftSubcarrier(i+cyclicPrefix_start,j);
    end
end

% Add cyclic prefix to the data
appendedCP = vertcat(cyclicPrefix_data, ifftSubcarrier);

% Convert to serial
dataOut = reshape(appendedCP,[numel(appendedCP),1]);

end
