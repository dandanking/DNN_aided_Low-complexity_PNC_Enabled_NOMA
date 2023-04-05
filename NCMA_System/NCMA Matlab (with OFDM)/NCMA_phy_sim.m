%% Use GNURadio's output to generate PHY map
%  Usage (fed correct input files;):
%  (1) Mode 1: from raw samples;
%  (2) Mode 2: from fft output;
%  (3) Mode 3: from demod output;
%
%  author: lzyou@inc.cuhk.edu.hk
%  date: June 1, 2013
%

clear all;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
global DECODER_TYPE
fft_size=64;
data_tones=48;
padding_bytes=1;
nsym=64;
nbits=nsym*data_tones;

Num_data_symbols=nsym;

DECODER_TYPE = 'soft';
mode = 2; %FIXME: change the code to be workable when mode=3
noise = 0.0000011;
% ------------------------------------------------------------------------%
%                       Self-Defined Part                                 %
% ------------------------------------------------------------------------%

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

%FIXME: normalized factor should be calculated automatically
% lzyou note (Nov 14,2013): according to the data, we should use 1.
% However, the performance when alpha=0.5 is bad. It is not consist with gnuradio results.
% Maybe gnuradio also use norm=8? Anyway, we use old setup, and will
% collect new data with new setup.
norm=8;
txA = read_complex_binary('ncma_txdata1/txA_64S.dat');
txA_TS_time = txA(161+16:160+80);
TS2_freq_a = fft(txA_TS_time);
TS2_freq_a = fftshift(TS2_freq_a)/norm;

txB = read_complex_binary('ncma_txdata1/txB_64S.dat');
txB_TS_time = txB(241+16:240+80);
TS2_freq_b = fftshift(fft(txB_TS_time))/norm;

null_encoded_bits=-10*ones(length(txX_encoded_bits),1);

%files={'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB'}; %,'0126_06_105dB'};
%files={'0126_06_105dB'};
%files={'0125_08_85dB'};

%files={'awgn_90_75dB'};
%files={'awgn_90_5dB','awgn_90_6dB','awgn_90_7dB','awgn_90_8dB','awgn_90_9dB'}; %,'awgn_90_95dB','awgn_90_10dB','awgn_90_105dB'};
%files={'awgn_cfo_flat_5dB','awgn_cfo_flat_6dB','awgn_cfo_flat_7dB','awgn_cfo_flat_8dB','awgn_cfo_flat_9dB'};
%files={'awgn_flat_5dB','awgn_flat_6dB','awgn_flat_7dB','awgn_flat_8dB','awgn_flat_9dB'};
%files={'awgn_45_5dB','awgn_45_6dB','awgn_45_7dB','awgn_45_8dB','awgn_45_9dB','awgn_0_5dB','awgn_0_6dB','awgn_0_7dB','awgn_0_8dB','awgn_0_9dB'};
%files={'75dB_75dB','75dB_85dB','75dB_95dB','75dB_105dB','75dB_115dB','75dB_125dB','75dB_135dB'};
files={'95dB_75dB','95dB_85dB','95dB_95dB','95dB_105dB','95dB_135dB'};
%files={'5dB','6dB','7dB','8dB','9dB'};
%files={'75dB_75dB','75dB_95dB','75dB_105dB','75dB_135dB'};

for fileIndex=1:length(files)

tag=files{fileIndex};

%ss=strcat('ncma_demod_data/awgn_flat_unbalanced_75dB/',tag);
%ss=strcat('ncma_demod_data/awgn_flat_unbalanced_95dB/',tag);
%ss=strcat('ncma_demod_data/awgn_random_degree/',tag);
%ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_75dB/',tag);
ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_95dB/',tag);
ss_in=strcat(ss,'/rx-ofdmdemod.dat');
ss_flag=strcat(ss,'/rx-ofdmdemod.datb');

fprintf(' Opening file: %s \n', tag);


% Mode 1: from raw samples
%FIXME: implement mode 1, and implement SNR calculation

%-------------------------------------------------------------------------%
%ss_in = '../pnc/rawofdm/examples/rx-ofdmdemod.dat';
%ss_flag = '../pnc/rawofdm/examples/rx-ofdmdemod.datb';

%ss_in = '/home/lzyou/pnc/rpnc/examples/logs/rx-fft.dat';
%ss_flag = '/home/lzyou/pnc/rpnc/examples/logs/rx-sampler.datb';
%-------------------------------------------------------------------------%


% Mode 2: begin from fft output
demod_in = read_complex_binary(ss_in);
demod_flag = read_char_binary(ss_flag);
pos_list=find(demod_flag>0);

num_pkts = length(pos_list);
hardRawBER     = zeros(num_pkts,3);
hardBER        = zeros(num_pkts,3);
hardPhyRawMap  = zeros(num_pkts,3);
phyRawMap      = zeros(num_pkts,3);
phyBridgeMap   = zeros(num_pkts,3);
phyTwoPhaseMap = zeros(num_pkts,3);

softBER = zeros(num_pkts,3);
softTwoPhaseBER = zeros(num_pkts,3);
snrMap = zeros(num_pkts,2);

A_INDEX=1; B_INDEX=2; X_INDEX=3;

% Mode 3: begin from soft demod output
% load GNURadio's soft output
%{
rxX=read_char_binary('../pnc-bits.dat');
rxA=read_char_binary('../pnc-a-bits.dat');
rxB=read_char_binary('../pnc-b-bits.dat');
%}



for jj=1:length(pos_list)
    
if mode == 2

    pos = pos_list(jj);
    rxdata = demod_in((pos-1)*fft_size+1:(pos+nsym-1+2)*fft_size);
    
    rxTS2A = rxdata(1:64);
    rxTS2B = rxdata(65:128);
    energyA = sum(abs(ifft(rxTS2A)).^2)/length(rxTS2A);
    energyB = sum(abs(ifft(rxTS2B)).^2)/length(rxTS2B);
    snrMap(jj,A_INDEX) = 10*log10(energyA/noise);
    snrMap(jj,B_INDEX) = 10*log10(energyB/noise);
    
    rxRawData = rxdata(129:end);
    
    H_a = rxTS2A ./ TS2_freq_a;
    H_b = rxTS2B ./ TS2_freq_b;
    
    NULL_INDEX = [1:6 33 60:64];
    H_a(NULL_INDEX) = 0;
    H_b(NULL_INDEX) = 0;
    
elseif mode == 3
    % Mode 3: begin from soft demod output
    % load GNURadio's soft output
    rawXBits = rxX((jj-1)*nbits+1:jj*nbits);
    rawABits = rxA((jj-1)*nbits+1:jj*nbits);
    rawBBits = rxB(((jj-1)*nbits+1:jj*nbits));

    rawHardXBits = rxX((jj-1)*nbits+1:jj*nbits)/255;
    rawHardABits = rxA((jj-1)*nbits+1:jj*nbits)/255;
    rawHardBBits = rxB(((jj-1)*nbits+1:jj*nbits))/255;
else
    fprintf(' Please specify the correct mode! \n');
end


% ------------------------------------------------------------------------%
HARD_INDEX=1;
SOFT_INDEX=2;
ITER_INDEX=3;

%FLAG = [1 1 0];
%FLAG = [0 1 0];
FLAG = [1 1 1];
%assert(FLAG(SOFT_INDEX)==FLAG(ITER_INDEX));

% ------------------------------------------------------------------------%
%                          Hard Decoder                                   %
% ------------------------------------------------------------------------%

if FLAG(HARD_INDEX)
    DECODER_TYPE = 'hard';
    knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
    RawHardBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
    rawHardXBits = RawHardBits(:,X_INDEX);
    rawHardABits = RawHardBits(:,A_INDEX);
    rawHardBBits = RawHardBits(:,B_INDEX);

    XBits = phy_viterbi_decoder(phy_de_interleaver(rawHardXBits),DECODER_TYPE);
    ABits = phy_viterbi_decoder(phy_de_interleaver(rawHardABits),DECODER_TYPE);
    BBits = phy_viterbi_decoder(phy_de_interleaver(rawHardBBits),DECODER_TYPE);
    % remove padding bits
    XBits=XBits(1:end-8);
    ABits=ABits(1:end-8);
    BBits=BBits(1:end-8);
    
    rawBerX = sum(abs(txX_encoded_bits-rawHardXBits))/length(rawHardXBits);
    rawBerA = sum(abs(txA_encoded_bits-rawHardABits))/length(rawHardABits);
    rawBerB = sum(abs(txB_encoded_bits-rawHardBBits))/length(rawHardBBits);
    hardRawBER(jj,[X_INDEX A_INDEX B_INDEX]) = [rawBerX rawBerA rawBerB];
    
    [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits,txX_source_bits,1);
    %fprintf(' rxCRCX=%s calCRCX=%s \n',data_bin2hex(rx_crc32_poly,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits,txA_source_bits,0);
    %fprintf(' %d: rxCRCA=%s calCRCA=%s ',jj,data_bin2hex(rx_crc32,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits,txB_source_bits,0);
    %fprintf(' rxCRCB=%s calCRCB=%s ',data_bin2hex(rx_crc32,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    fprintf('  hard %d: okX=%d berX=[%f %f] | okA=%d berA=[%f %f] | okB=%d berB=[%f %f] | snrA=%f snrB=%f \n', jj, okX, berX, rawBerX, okA, berA, rawBerA, okB, berB, rawBerB, snrMap(jj,A_INDEX), snrMap(jj,B_INDEX));
    
    hardPhyRawMap(jj,[X_INDEX A_INDEX B_INDEX]) = [okX okA okB];
    hardBER(jj,[X_INDEX A_INDEX B_INDEX]) = [berX berA berB];

end

% ------------------------------------------------------------------------%

% ------------------------------------------------------------------------%
%                          Soft Decoder                                   %
% ------------------------------------------------------------------------%
%%{

if FLAG(SOFT_INDEX)
    DECODER_TYPE = 'soft';
    knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
    RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);    
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
    
    [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits,txX_source_bits,1);
    %fprintf(' rxCRCX=%s calCRCX=%s \n',data_bin2hex(rx_crc32_poly,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits,txA_source_bits,0);
    %fprintf(' %d: rxCRCA=%s calCRCA=%s ',jj,data_bin2hex(rx_crc32,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits,txB_source_bits,0);
    %fprintf(' rxCRCB=%s calCRCB=%s ',data_bin2hex(rx_crc32,'left-msb'),data_bin2hex(cal_crc32,'left-msb'));
    fprintf('  soft %d: okX=%d berX=%f | okA=%d berA=%f | okB=%d berB=%f | snrA=%f snrB=%f  \n', jj, okX, berX, okA, berA, okB, berB, snrMap(jj,A_INDEX), snrMap(jj,B_INDEX));
    
    phyRawMap(jj,[X_INDEX A_INDEX B_INDEX]) = [okX okA okB];
    phyBridgeMap(jj,[X_INDEX A_INDEX B_INDEX]) = [okX okA okB];
    softBER(jj,[X_INDEX A_INDEX B_INDEX]) = [berX berA berB];
    
    %%{
    %-------------------------------------------------------------------------%
    if (okX&&okA&&~okB) || (okX&&okB&&~okA)
        berX = 0; berA = 0; berB = 0;
        okX = 1; okA = 1; okB = 1;
        phyBridgeMap(jj,[X_INDEX A_INDEX B_INDEX]) = [okX okA okB];
    else %-two-pass decoding
        if FLAG(ITER_INDEX)
            knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
            if okX&&(~okA)&&(~okB)
                knownTxRawData(:,X_INDEX) = txX_encoded_bits;
                RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawABits = RawBits(:,A_INDEX); rawBBits = RawBits(:,B_INDEX);
                ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
                BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
                [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
                [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
                fprintf('  soft+%d: okX=%d berX=%f | okA=%d berA=%f | okB=%d berB=%f  \n', jj, okX, berX, okA, berA, okB, berB);
                if okA || okB
                    okA = 1; okB = 1;
                    berA = 0; berB = 0;
                    %fprintf(' =%d: okA=%d okB=%d okX=%d \n', jj, okA, okB, 1);
                end
            elseif okA&&(~okX)&&(~okB)
                knownTxRawData(:,A_INDEX) = txA_encoded_bits;
                RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawXBits = RawBits(:,X_INDEX); rawBBits = RawBits(:,B_INDEX);
                XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),DECODER_TYPE);
                BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
                [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits(1:end-8),txX_source_bits,1);
                [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
                fprintf('  soft+%d: okX=%d berX=%f | okA=%d berA=%f | okB=%d berB=%f  \n', jj, okX, berX, okA, berA, okB, berB);
                if okX || okB
                    %fprintf(' =%d: okA=%d okB=%d okX=%d \n', jj, 1, okB, okX);
                    okX = 1; okB = 1;
                    berX = 0; berB = 0;
                end
            elseif okB&&(~okX)&&(~okA)
                knownTxRawData(:,B_INDEX) = txB_encoded_bits;
                RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawXBits = RawBits(:,X_INDEX); rawABits = RawBits(:,A_INDEX);
                XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),DECODER_TYPE);
                ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
                [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits(1:end-8),txX_source_bits,1);
                [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
                fprintf('  soft+%d: okX=%d berX=%f | okA=%d berA=%f | okB=%d berB=%f  \n', jj, okX, berX, okA, berA, okB, berB);
                if okX || okA
                    %fprintf(' =%d: okA=%d okB=%d okX=%d \n', jj, okA, 1, okX);
                    okX = 1; okA = 1;
                    berX = 0; berA = 0;
                end
                
            end
        end % end of if two-phase decoding
        
    end % end of PHY bridging

    phyTwoPhaseMap(jj,[X_INDEX A_INDEX B_INDEX]) = [okX okA okB];
    softTwoPhaseBER(jj,[X_INDEX A_INDEX B_INDEX]) = [berX berA berB];
end %end of if soft

end


xright=length(find(hardPhyRawMap(:,X_INDEX)));
aright=length(find(hardPhyRawMap(:,A_INDEX)));
bright=length(find(hardPhyRawMap(:,B_INDEX)));
fprintf(' Hard: %d %d %d \n',xright,aright,bright);

xright=length(find(phyRawMap(:,X_INDEX)));
aright=length(find(phyRawMap(:,A_INDEX)));
bright=length(find(phyRawMap(:,B_INDEX)));
fprintf(' Soft-Raw: %d %d %d \n',xright,aright,bright);

xright=length(find(phyBridgeMap(:,X_INDEX)));
aright=length(find(phyBridgeMap(:,A_INDEX)));
bright=length(find(phyBridgeMap(:,B_INDEX)));
fprintf(' Soft-Bridge: %d %d %d \n',xright,aright,bright);

xright=length(find(phyTwoPhaseMap(:,X_INDEX)));
aright=length(find(phyTwoPhaseMap(:,A_INDEX)));
bright=length(find(phyTwoPhaseMap(:,B_INDEX)));
fprintf(' Soft-TwoPhase: %d %d %d \n',xright,aright,bright);

clear demod_in;
clear demod_flag;
ss=strcat('matlab_data/',tag);
ss=strcat(ss,'_maps.mat');
%save(ss,'hardPhyRawMap','phyRawMap','phyBridgeMap','phyTwoPhaseMap');
save(ss);


end
