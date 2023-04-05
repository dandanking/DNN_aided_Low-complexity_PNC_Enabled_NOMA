% Long Training Symbol channel estimation function
% by lulu 20100821
% by HB 20100821

% -target- %
% To use the Long training symbol to decode the channel matrix inv(H)
% order-changned data FFT symbol.

% modified by Taotao 2010-11-25

function [Delta_Theta] = My_CFO_Est (Train_Symbol)

N = 64; % the number of points for FFT
Index = 8; % we choose the k-th index for the FFT symbol begining.
% long train symbol input, to copy only the long training symbol for
% processing.
Train_Time = Train_Symbol(1:80);

Train_Time_Pt1 = Train_Time(Index+1 : Index+N/2); % first half of 32 samples
Train_Time_Pt2 = Train_Time(Index+N/2+1 : Index+N); % second half of 32 samples

Pt1_Theta = atan2( imag(Train_Time_Pt1), real(Train_Time_Pt1) ); % the phase of the first partition 128x1
Pt2_Theta = atan2( imag(Train_Time_Pt2), real(Train_Time_Pt2) ); % the phase of the first partition 128x1

% Delta_Theta = mean(Pt2_Theta - Pt1_Theta) / 128; % only one value. This is the averaged Delta_Theta per sample
% -Delta_Theta Estimation finished- %
Delta_Theta32 = (Pt2_Theta - Pt1_Theta) / 32;
Delta_Theta_Sort = sort(Delta_Theta32);
Delta_Theta = Delta_Theta_Sort(17); % we use the middle number to replace the averaged number


%% new LTS CFO estimation --2011-12-05 ----added by lulu
%% new LTS CFO estimation --2012-03-08 ----modified by lzyou
%------------------------------------------------------%
%     long_product = conj(Train_Time_Pt1) .* Train_Time_Pt2;
%     long_sum = sum(long_product);
%     phase = atan2(imag(long_sum), real(long_sum))/32;
%     Delta_Theta1 =  phase;
%------------------------------------------------------%
%% -end STS CFO estimation
end