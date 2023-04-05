%% Use GNURadio's output to generate SIC map
%
%  Input:
%    * signal after fft
%
%  Output:
%    * packet decoding statistics
%
%  author: lzyou@inc.cuhk.edu.hk
%  date: June 1, 2013
%

clear all;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
global DECODER_TYPE
global alpha
fft_size=64;
data_tones=48;
padding_bytes=1;
nsym=64;
nbits=nsym*data_tones;

Num_data_symbols=nsym;

%DECODER_TYPE = 'hard';
DECODER_TYPE = 'soft';
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

null_encoded_bits=-10*ones(length(txX_encoded_bits),1);

%FIXME: normalized factor should be calculated automatically
norm=8;
txA = read_complex_binary('ncma_txdata1/txA_64S.dat');
txA_TS_time = txA(161+16:160+80);
TS2_freq_a = fft(txA_TS_time);
TS2_freq_a = fftshift(TS2_freq_a)/norm;

txB = read_complex_binary('ncma_txdata1/txB_64S.dat');
txB_TS_time = txB(241+16:240+80);
TS2_freq_b = fftshift(fft(txB_TS_time))/norm;


files={'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
files={'0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB'};
files={'0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB'};
files={'0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB'};
files={'0126_06_105dB'};
files={'75dB_75dB','75dB_85dB','75dB_95dB','75dB_105dB','75dB_115dB','75dB_125dB','75dB_135dB'};
files={'95dB_75dB','95dB_85dB','95dB_95dB','95dB_105dB','95dB_135dB'};
%files={'5dB','6dB','7dB','8dB','9dB'};
%files={'75dB_75dB','75dB_95dB','75dB_105dB','75dB_135dB'};

%alpha_list = {0.025, 0.05, 0.1, 0.25};
alpha_list = {0.05};
%alpha_list = {0.1, 0.25, 0.5};

for fileIndex=1:length(files)

tag=files{fileIndex};

%ss=strcat('ncma_demod_data/awgn_flat_unbalanced_75dB/',tag);
%ss=strcat('ncma_demod_data/awgn_flat_unbalanced_95dB/',tag);
%ss=strcat('ncma_demod_data/awgn_random_degree/',tag);
%ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_75dB/',tag);
ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_95dB/',tag);
ss_in=strcat(ss,'/rx-ofdmdemod.dat');
ss_flag=strcat(ss,'/rx-ofdmdemod.datb');

fprintf(' \n Opening file: %s \n', tag);

DEBUG = 1;

%-------------------------------------------------------------------------%
%ss_in = '../pnc/rawofdm/examples/rx-ofdmdemod.dat';
%ss_flag = '../pnc/rawofdm/examples/rx-ofdmdemod.datb';

%ss_in = '/home/lzyou/pnc/rpnc/examples/logs/rx-fft.dat';
%ss_flag = '/home/lzyou/pnc/rpnc/examples/logs/rx-sampler.datb';
%-------------------------------------------------------------------------%

for alphaIndex=1:length(alpha_list)
    alpha = alpha_list{alphaIndex};
    
    %-----------------------------------------%
    %       Specify HARD or SOFT here         %
    %-----------------------------------------%
    A_INDEX=1; B_INDEX=2; X_INDEX=3;
    HARD_INDEX=1; SOFT_INDEX=2;
    STR_INDEX = {'hard', 'soft'};
    INDEX = [HARD_INDEX, SOFT_INDEX];
    VALUE = [0, 1];
    FLAG = zeros(length(INDEX),1);
    FLAG(INDEX) = VALUE;
    %-----------------------------------------%

for ff=1:length(INDEX)
    
    demod_in = read_complex_binary(ss_in);
    demod_flag = read_char_binary(ss_flag);
    pos_list=find(demod_flag>0);
    
    num_pkts = length(pos_list);
    hardRawBER     = zeros(num_pkts,3);
    hardBER        = zeros(num_pkts,3);
    hardSIC1BER    = zeros(num_pkts,3);
    hardSIC2BER    = zeros(num_pkts,3);
    softSIC1BER    = zeros(num_pkts,3);
    softSIC2BER    = zeros(num_pkts,3);
    hardPhyRawMap  = zeros(num_pkts,3);
    hardPhySIC1Map = zeros(num_pkts,3);
    hardPhySIC2Map = zeros(num_pkts,3);
    softPhyRawMap  = zeros(num_pkts,3);
    softPhySIC1Map = zeros(num_pkts,3);
    softPhySIC2Map = zeros(num_pkts,3);
    
    softBER = zeros(num_pkts,3);
    snrMap = zeros(num_pkts,2);
    multiUser = zeros(num_pkts,1); % indicator of two-user access

for jj=1:length(pos_list)
    
    pos = pos_list(jj);
    rxdata = demod_in((pos-1)*fft_size+1:(pos+nsym-1+2)*fft_size);
    
    rxTS2A = rxdata(1:64);
    rxTS2B = rxdata(65:128);
    energyA = sum(abs(ifft(rxTS2A)).^2)/length(rxTS2A);
    energyB = sum(abs(ifft(rxTS2B)).^2)/length(rxTS2B);
    snrMap(jj,A_INDEX) = 10*log10(energyA/noise);
    snrMap(jj,B_INDEX) = 10*log10(energyB/noise);
    if snrMap(jj,A_INDEX) > 1 && snrMap(jj,B_INDEX) > 1
        multiUser(jj) = 1;
    end
    
    rxRawData = rxdata(129:end);
    
    H_a = rxTS2A ./ TS2_freq_a;
    H_b = rxTS2B ./ TS2_freq_b;
    
    NULL_INDEX = [1:6 33 60:64];
    H_a(NULL_INDEX) = 0;
    H_b(NULL_INDEX) = 0;

    if multiUser(jj) == 0
        if DEBUG
            fprintf('  Single User: SNR_A=%f SNR_B=%f \n',snrMap(jj,A_INDEX),snrMap(jj,B_INDEX));
        end
        continue;
    end

    if ff == HARD_INDEX
        DECODER_TYPE = 'hard';
    elseif ff == SOFT_INDEX
        DECODER_TYPE = 'soft';
    end
    
    if FLAG(ff)
        knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
        [RawBits,~,~] = SIC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
        rawABits = RawBits(:,A_INDEX);
        rawBBits = RawBits(:,B_INDEX);
        
        ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
        BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
        % remove padding bits
        ABits=ABits(1:end-8);
        BBits=BBits(1:end-8);
        
        [okA,berA,~,~] = mac_crc32_wrapper(ABits,txA_source_bits,0);
        [okB,berB,~,~] = mac_crc32_wrapper(BBits,txB_source_bits,0);
        
        if DEBUG
            fprintf('  %s %d: okA=%d berA=%f | okB=%d berB=%f | snrA=%f snrB=%f  \n', STR_INDEX{ff}, jj, okA, berA, okB, berB, snrMap(jj,A_INDEX), snrMap(jj,B_INDEX));
        end
        
        if ff == HARD_INDEX
            rawBerA = sum(abs(txA_encoded_bits-rawABits))/length(rawABits);
            rawBerB = sum(abs(txB_encoded_bits-rawBBits))/length(rawBBits);
            hardRawBER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 rawBerA rawBerB];
            
            hardPhyRawMap(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            hardBER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        elseif ff == SOFT_INDEX
            softPhyRawMap(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            softBER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        end
        
        sic_okA = okA; sic_okB = okB;
        % SIC1 Algorithm: strong user first
        if okA || okB
            if okA && snrMap(jj,A_INDEX) < snrMap(jj,B_INDEX)
                okA = 0;
            elseif okB && snrMap(jj,A_INDEX) > snrMap(jj,B_INDEX)
                okB = 0;
            elseif okA && snrMap(jj,A_INDEX) > snrMap(jj,B_INDEX)
                knownTxRawData(:,A_INDEX) = txA_encoded_bits;
                [RawBits,~,SIC_B_debug] = SIC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawBBits = RawBits(:,B_INDEX);
                BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
                [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
                if DEBUG
                    fprintf('  %s-SIC1-%d: okA=%d berA=%f | okB=%d berB=%f  \n', STR_INDEX{ff}, jj, okA, berA, okB, berB);
                end
            elseif okB && snrMap(jj,A_INDEX) < snrMap(jj,B_INDEX)
                knownTxRawData(:,B_INDEX) = txB_encoded_bits;
                [RawBits,~,SIC_A_debug] = SIC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawABits = RawBits(:,A_INDEX);
                ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
                [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
                if DEBUG
                    fprintf('  %s-SIC1-%d: okA=%d berA=%f | okB=%d berB=%f  \n', STR_INDEX{ff}, jj, okA, berA, okB, berB);
                end
            end
        end
        
        if ff == HARD_INDEX
            hardPhySIC1Map(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            hardSIC1BER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        elseif ff == SOFT_INDEX
            softPhySIC1Map(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            softSIC1BER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        end
        
        % SIC2 Algorithm: parallel SIC
        okA = sic_okA; okB = sic_okB;
        knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
        if okA&&(~okB)
            knownTxRawData(:,A_INDEX) = txA_encoded_bits;
            [RawBits,~,SIC_B_debug] = SIC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
            rawBBits = RawBits(:,B_INDEX);
            BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
            [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
            if DEBUG
                fprintf('  %s-SIC2-%d: okA=%d berA=%f | okB=%d berB=%f  \n', STR_INDEX{ff}, jj, okA, berA, okB, berB);
            end
        elseif okB&&(~okA)
            knownTxRawData(:,B_INDEX) = txB_encoded_bits;
            [RawBits,~,SIC_A_debug] = SIC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
            rawABits = RawBits(:,A_INDEX);
            ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
            [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
            if DEBUG
                fprintf('  %s-SIC2-%d: okA=%d berA=%f | okB=%d berB=%f  \n', STR_INDEX{ff}, jj, okA, berA, okB, berB);
            end
        end
        
        if ff == HARD_INDEX
            hardPhySIC2Map(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            hardSIC2BER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        elseif ff == SOFT_INDEX
            softPhySIC2Map(jj,[X_INDEX A_INDEX B_INDEX]) = [0 okA okB];
            softSIC2BER(jj,[X_INDEX A_INDEX B_INDEX]) = [1 berA berB];
        end
    end 
end % end of pos_list

if FLAG(ff) && ff == HARD_INDEX
    fprintf(' alpha = %f \n', alpha);
    
    aright=length(find(hardPhyRawMap(:,A_INDEX)));
    bright=length(find(hardPhyRawMap(:,B_INDEX)));
    fprintf(' Hard-Raw: %d %d \n',aright,bright);
    
    aright=length(find(hardPhySIC1Map(:,A_INDEX)));
    bright=length(find(hardPhySIC1Map(:,B_INDEX)));
    fprintf(' Hard-SIC1: %d %d \n',aright,bright);
    
    aright=length(find(hardPhySIC2Map(:,A_INDEX)));
    bright=length(find(hardPhySIC2Map(:,B_INDEX)));
    fprintf(' Hard-SIC2: %d %d \n',aright,bright);
elseif FLAG(ff) && ff == SOFT_INDEX
    fprintf(' alpha = %f \n', alpha);
    
    aright=length(find(softPhyRawMap(:,A_INDEX)));
    bright=length(find(softPhyRawMap(:,B_INDEX)));
    fprintf(' Soft-Raw: %d %d \n',aright,bright);
    
    aright=length(find(softPhySIC1Map(:,A_INDEX)));
    bright=length(find(softPhySIC1Map(:,B_INDEX)));
    fprintf(' Soft-SIC1: %d %d \n',aright,bright);
    
    aright=length(find(softPhySIC2Map(:,A_INDEX)));
    bright=length(find(softPhySIC2Map(:,B_INDEX)));
    fprintf(' Soft-SIC2: %d %d \n',aright,bright);
end

if FLAG(ff)
    clear demod_in;
    clear demod_flag;
    ss=strcat('matlab_data/',tag);
    ss=strcat(ss,'_');
    ss=strcat(ss,STR_INDEX{ff});
    ss=strcat(ss,'SIC_alpha_');
    ss=strcat(ss,num2str(alpha));
    ss=strcat(ss,'_map.mat');
    save(ss);
end

end % end of INDEX

end % end of alpha

end