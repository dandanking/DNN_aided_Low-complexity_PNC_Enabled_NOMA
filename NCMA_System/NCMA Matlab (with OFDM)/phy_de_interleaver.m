% de_interleaver for 802.11 receiver 
% TODO: implement the second permutation for 16QAM and 64QAM
% 2014-7-29 

function data_deinterleaved = phy_de_interleaver(data_interleaved)
global MOD
if strcmp(MOD,'BPSK')
    cbps = 48;
elseif strcmp(MOD,'QPSK')
    cbps = 48*2;
%elseif strcmp(MOD,'16QAM')
%    cbps = 48*4;
%elseif strcmp(MOD,'64QAM')
%    cbps = 48*6;
end

din = reshape (data_interleaved, cbps, []);
for i = 0 : cbps-1
    k = cbps/16*mod(i,16) + floor(i/16);
    dtemp(i+1,:) = din(k+1,:);
end

data_deinterleaved = reshape(dtemp, 1, []);
data_deinterleaved = data_deinterleaved';


