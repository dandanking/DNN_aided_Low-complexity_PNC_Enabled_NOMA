
close all;
clear all;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
global DECODER_TYPE
fft_size=64;
data_tones=48;
padding_bytes=1;

nsym = 16;
Num_data_symbols=nsym;
nbits=nsym*data_tones;
RATE=1;
global MOD
if RATE <= 2
  MOD = 'BPSK';
elseif RATE <= 4
  MOD = 'QPSK';
end
DECODER_TYPE = 'soft';

% ------------------------------------------------------------------------%
%                          Tx Data Part                                   %
% ------------------------------------------------------------------------%
strA = sprintf('ncma_txdata2/txA_BPSK_%dS_intrlv.datb',Num_data_symbols);
strB = sprintf('ncma_txdata2/txB_BPSK_%dS_intrlv.datb',Num_data_symbols);
txA_encoded_bits=read_char_binary(strA);
txB_encoded_bits=read_char_binary(strB);
txX_encoded_bits=bitxor(txA_encoded_bits,txB_encoded_bits);
txA_source_bits=phy_viterbi_decoder(phy_de_interleaver(txA_encoded_bits),'hard');
txB_source_bits=phy_viterbi_decoder(phy_de_interleaver(txB_encoded_bits),'hard');
txX_source_bits=phy_viterbi_decoder(phy_de_interleaver(txX_encoded_bits),'hard');
% remove padding bits
txA_source_bits=txA_source_bits(1:end-8);
txB_source_bits=txB_source_bits(1:end-8);
txX_source_bits=txX_source_bits(1:end-8);

strA = sprintf('ncma_txdata2/txA_BPSK_%dS.dat',Num_data_symbols);
strB = sprintf('ncma_txdata2/txB_BPSK_%dS.dat',Num_data_symbols);

norm=8;
txA = read_complex_binary(strA);
txA_TS_time = txA(161+16:160+80);
TS2_freq_a = fft(txA_TS_time);
TS2_freq_a = fftshift(TS2_freq_a)/norm;

txB = read_complex_binary(strB);
txB_TS_time = txB(241+16:240+80);
TS2_freq_b = fftshift(fft(txB_TS_time))/norm;

% Mode 1: begin after sync
%ss_in='debug_data/rx-sigmix.dat';
%ss_flag='debug_data/rx-sync.datb';

% Mode 2: begin after sampler
ss_in=('../matlab-mimo/new_data_mac/rx-fft-bpsk.dat');
ss_flag=('../matlab-mimo/new_data_mac/rx-sampler-bpsk.datb');

sampler_out = read_complex_binary(ss_in);
sampler_flag = read_char_binary(ss_flag);
pos_list=find(sampler_flag>0);

A_INDEX=1; B_INDEX=2; X_INDEX=3;
null_encoded_bits=-10*ones(length(txX_encoded_bits),1);

for jj = 1:1%length(pos_list)
    pos = pos_list(jj);
    rxdata = sampler_out((pos-1)*fft_size+1:(pos+nsym-1+2)*fft_size);
    
    rxTS2A = rxdata(1:64);
    rxTS2B = rxdata(65:128);
    rxRawData = rxdata(129:end);
    
    H_a = rxTS2A ./ TS2_freq_a;
    H_b = rxTS2B ./ TS2_freq_b;
    
    NULL_INDEX = [1:6 33 60:64];
    H_a(NULL_INDEX) = 0;
    H_b(NULL_INDEX) = 0;
    
    knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
    [RawBits,debug] = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);    
    rawXBits = RawBits(:,X_INDEX);
    rawABits = RawBits(:,A_INDEX);
    rawBBits = RawBits(:,B_INDEX);

    XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),DECODER_TYPE);
    ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
    BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
    % remove padding bits
    XBits=XBits(1:end-8);
    ABits=ABits(1:end-8);
    BBits=BBits(1:end-8);
    
    fprintf('[%d %d %d] \n',sum(abs(txA_source_bits-ABits)),sum(abs(txB_source_bits-BBits)),sum(abs(txX_source_bits-XBits)));

end