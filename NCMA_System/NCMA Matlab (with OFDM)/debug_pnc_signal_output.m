
close all;
clear all;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
global DECODER_TYPE
symbol_size=80;
fft_size=64;
data_tones=48;
padding_bytes=1;
nsym=64;
nbits=nsym*data_tones;

Num_data_symbols=nsym;

DECODER_TYPE = 'hard';

RATE = 1;
global MOD
if RATE <= 2
  MOD = 'BPSK';
elseif RATE <= 4
  MOD = 'QPSK';
end

% Mode: 
%   0: raw samples
%   1: raw samples with indicator (without cp cut & fft)
mode = 0; %[0,1]

% ------------------------------------------------------------------------%
%                          Tx Data Part                                   %
% ------------------------------------------------------------------------%
txA_encoded_bits=read_char_binary('ncma_txdata1/txA_64S_intrlv.datb');
txB_encoded_bits=read_char_binary('ncma_txdata1/txB_64S_intrlv.datb');
txX_encoded_bits=bitxor(txA_encoded_bits,txB_encoded_bits);
txA_source_bits=phy_viterbi_decoder(phy_de_interleaver(txA_encoded_bits),'hard');
txB_source_bits=phy_viterbi_decoder(phy_de_interleaver(txB_encoded_bits),'hard');
txX_source_bits=phy_viterbi_decoder(phy_de_interleaver(txX_encoded_bits),'hard');
% remove padding bits
txA_source_bits=txA_source_bits(1:end-8);
txB_source_bits=txB_source_bits(1:end-8);
txX_source_bits=txX_source_bits(1:end-8);

norm=8;
txA = read_complex_binary('ncma_txdata1/txA_64S.dat');
txA_TS_time = txA(161+16:160+80);
TS2_freq_a = fftshift(fft(txA_TS_time))/norm;

txB = read_complex_binary('ncma_txdata1/txB_64S.dat');
txB_TS_time = txB(241+16:240+80);
TS2_freq_b = fftshift(fft(txB_TS_time))/norm;

% Mode 0: raw samples
if mode == 0
    ss_in = 'debug_data/ncma_samples_0126_06/ncma_samples_0126_06_3pkts.dat';
    rxsig = read_complex_binary(ss_in);
end

% Mode 1: begin after sync
if mode == 1
    ss_in='debug_data/rx-sigmix.dat';
    ss_flag='debug_data/rx-sync.datb';

    sampler_out = read_complex_binary(ss_in);
    sampler_flag = read_char_binary(ss_flag);
    pos_list=find(sampler_flag>0);
end

A_INDEX=1; B_INDEX=2; X_INDEX=3;
null_encoded_bits=-10*ones(length(txX_encoded_bits),1);

cp_pos = -1;
for jj = 1:1
    
    if mode == 0
        % debug % TODO: implement frame sync function
        % pkt1: 84273; pkt2: 384427; pkt3: 684419
        pos = 84273;
        cp_pos = 6;
        RawData = rxsig(pos+1+cp_pos:pos+cp_pos+(nsym+2)*symbol_size);
    end
    
    if mode == 1
        pos = pos_list(jj);
        RawData = sampler_out(pos+1:pos+(nsym+2)*symbol_size);
    end
    
    if mode == 0 || mode == 1
        rxdata = zeros(64/80*length(RawData),1);
        m = 1;
        for i = 1 : 80 : length(RawData)
            timeRawData = RawData(i:i+fft_size-1);  % cp is already cut
            freqRawData = fft(timeRawData);
            freqRawData = fftshift(freqRawData);
            rxdata((m-1)*fft_size+1 : m*fft_size) = freqRawData;
            m = m + 1;
        end
    end
        
    
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
   
    errXBits = sum(abs(txX_source_bits-XBits));
    errABits = sum(abs(txA_source_bits-ABits));
    errBBits = sum(abs(txB_source_bits-BBits));
    fprintf('[ mode=%d cp_pos=%d file=%s] \n',mode,cp_pos,ss_in);
    fprintf('Error Bits (X,A,B) = [%d %d %d] [%f %f %f] \n',errXBits,errABits,errBBits,errXBits/length(XBits),errABits/length(ABits),errBBits/length(BBits));

end