% Long Training Symbol channel estimation function

% -target- %
% To use the Long training symbol to decode the channel matrix inv(H)
% order-changned data FFT symbol.

% modified by Taotao 2010-11-25
% revised by lzyou at Nov 11, 2013

function [H] = phy_channel_estimation (Train_Symbol, Tx_Train_Symbol, CutPos)

N = 64; % the number of points for FFT
Index = CutPos; % we choose the k-th index for the FFT symbol begining.
% long train symbol input, to copy only the long training symbol for
% processing.

Train_Freq = fftshift(fft(Train_Symbol(Index+1 : Index+N)));
Tx_Train_Freq = fftshift(fft(Tx_Train_Symbol(16 + 1 : 16 + 64)));
H = Train_Freq./Tx_Train_Freq;

range = [1:6 33 60:64];
H(range) = 0;

end