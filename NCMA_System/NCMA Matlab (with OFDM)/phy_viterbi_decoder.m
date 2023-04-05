% Viterbi decoder for convolutional code
% 2011-1-17

% by lulu 23/11/2012 to add the soft decoding module
% edited by lulu and lzyou 19/3/2013: to change the padding zeros, and
% trackback depth (larger tblen is better but slower)

function output_symbol = phy_viterbi_decoder(input_symbol,type)
% CodeGen = [133 171]; %generator polynomials used for FTW OFDM
CodeGen = [117 155];   %used for Raw OFDM
K = 7;%constraint length 
Trellis = poly2trellis(K,CodeGen); 

tblen = 48; %traceback depth   
signal_encoded = [input_symbol; zeros(tblen*2,1);]; %pading zeros at the end 

if strcmp(type,'hard') == 1
    [signal_decode m p in] = vitdec(signal_encoded, Trellis, tblen, 'cont', 'hard');% viterbi hard decoding 
elseif strcmp(type, 'soft') == 1
    [signal_decode m p in] = vitdec(signal_encoded, Trellis, tblen, 'cont', 'soft', 8);% viterbi soft decoding to decode the symbols to 0-255 level, 255 = 2^8 - 1 
elseif strcmp(type, 'exact') == 1
    [signal_decode m p in] = vitdec(signal_encoded, Trellis, tblen, 'cont', 'soft', 8);% viterbi soft decoding to decode the symbols to 0-255 level, 255 = 2^8 - 1
else
    fprintf(' Viterbi Error: please specify decoder type \n');
end

output_symbol = signal_decode (tblen+1:tblen+length(input_symbol)/2);%removing zeros form the begining

