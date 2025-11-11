function [dataOut] = demodOFDM(dataIn,cpLen,ofdmSym)

%--------------------------------------------------------------------------
%
%              Demodulates OFDM symbols into serial data
%
%--------------------------------------------------------------------------
% Definition of input arguments 
%
% dataIn                         Input data vector
% cpLen                          Length of the cyclic prefix
% ofdmSym                        No. of ofdm symbols per subframe
%
%--------------------------------------------------------------------------
% Function returns: 
%
% dataOut                        Output time domain symbols
%
%--------------------------------------------------------------------------

% OFDM receiever reshapes serial data to parallel
parallelRx = reshape(dataIn, numel(dataIn)/ofdmSym, ofdmSym);
% Removing the cyclic Prefix
parallelRx(1:(cpLen), :) = [];

% Perform FFT
dataOut =  fft(parallelRx,[],2); 


end