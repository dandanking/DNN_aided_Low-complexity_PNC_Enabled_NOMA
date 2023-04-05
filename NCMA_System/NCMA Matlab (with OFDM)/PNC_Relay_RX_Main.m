%% lzyou: Oct. 13, 2014
%% Note: set RATE, MODE, DECODER_TYPE, index
%%     RATE: modulation + coding rate (use 8 rates as in 802.11)
%%     MODE: AWGN or REAL
%%     DECODER_TYPE: hard or soft
%%     Num_data_symbols: number of data symbols
%%     index: the beginning index of the frame (before STS)
%% TODO:
%%     develop a peak detector to avoid setting index manualy

%clear all;

Num_preamble_symbols = 4;
Num_data_symbols  =  16;  %[16,64,256]
bandwidth = 5; 
RATE = 1;

global MOD
if RATE <= 2
  MOD = 'BPSK';
elseif RATE <= 4
  MOD = 'QPSK';
end

global DECODER_TYPE
%DECODER_TYPE='soft';
DECODER_TYPE='hard';

strA = sprintf('ncma_txdata2/txA_%s_%dS_intrlv.datb',MOD,Num_data_symbols);
strB = sprintf('ncma_txdata2/txB_%s_%dS_intrlv.datb',MOD,Num_data_symbols);
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

strA = sprintf('ncma_txdata2/txA_%s_%dS.dat',MOD,Num_data_symbols);
strB = sprintf('ncma_txdata2/txB_%s_%dS.dat',MOD,Num_data_symbols);
txa = read_complex_binary(strA);
txb = read_complex_binary(strB);

%str = sprintf('benchmark_rxdata/awgn_random_phase_%s_AB_%dS_40db_1000Pkts.dat',MOD,Num_data_symbols);
str = sprintf('benchmark_rxdata/awgn_random_phase_%s_AB_%dS.dat',MOD,Num_data_symbols);
%rx = read_complex_binary(str);
MODE = 'AWGN';

str = sprintf('benchmark_rxdata/realchannel_%dM_%s_AB_%dS.dat',bandwidth,MOD,Num_data_symbols);
%str = sprintf('benchmark_rxdata/rx-sampler-logging.dat');
%rx = read_complex_binary(str);
rx = read_complex_binary('../matlab-mimo/new_data_mac/rx-sampler-logging-bpsk.dat');
MODE = 'REAL';

if strcmp(MODE,'AWGN')
    GAP = 3*length(txa);
elseif strcmp(MODE,'REAL')
    GAP= ((Num_data_symbols+Num_preamble_symbols)*2)*80;
    if length(rx) > 1000*80+1
      rx = rx(1000*80+1:end);
    else
      rx = rx;
    end
end

Num_frame = 10; %100;

berX_all=[];
berA_all=[];
berB_all=[];
measureEVMA=[];
measureEVMB=[];

for ff = 1 : Num_frame
    
    % ------------------------------------------------------------------- %
    %FIXME: develop rawofdm-compatible peak detector function
    %[user, index] = Find_Sync_Index(rx, 0);
    %index = 1064-160;
    %index=160-160;
    if Num_data_symbols == 16
        if strcmp(MODE,'AWGN')
            index = 1760-160;
        elseif strcmp(MODE,'REAL')
            index = 0;%4677-160;
        end
    elseif Num_data_symbols == 64
        if strcmp(MODE,'AWGN')
            index=5600-160;
        elseif strcmp(MODE,'REAL')
            index = 934-160;
        end
    elseif Num_data_symbols == 256
        if strcmp(MODE,'AWGN')
            index = 20960-160;
        elseif strcmp(MODE,'REAL')
            index = 4886-160;
        end
    else
        index = 1000*80;
    end
    % ------------------------------------------------------------------- %

    frame = rx(index+1:index+length(txa));
    
    rx = rx(GAP+1:end);

    [XBits, ABits, BBits,enXBits, enABits, enBBits,measureA,measureB] = PNC_Packet_Decoder(frame, Num_data_symbols, txa, txb); %rx

    [okX,berX] = mac_crc32_wrapper(XBits,txX_source_bits,1);
    [okA,berA] = mac_crc32_wrapper(ABits,txA_source_bits,0);
    [okB,berB] = mac_crc32_wrapper(BBits,txB_source_bits,0);
    
    berX_all=[berX_all,berX];
    berA_all=[berA_all,berA];
    berB_all=[berB_all,berB];
    measureEVMA=[measureEVMA, measureA];
    measureEVMB=[measureEVMB, measureB];
    
    fprintf(' ---- MODE=%s NSYM=%d DECODER_TYPE=%s \n', MODE, Num_data_symbols, DECODER_TYPE);
    if okX
        fprintf(' ff=%d XOR CRC32 OK, ', ff);
    else
        fprintf(' ff=%d XOR CRC32 Fail, ', ff);
    end
    
    if okA
        fprintf('A CRC32 OK, ');
    else
        fprintf('A CRC32 Fail, ');
    end

    if okB
        fprintf('B CRC32 OK\n');
    else
        fprintf('B CRC32 Fail\n');
    end
    
    fprintf(' ff=%d berX=%f berA=%f berB=%f \n', ff, berX, berA, berB);

end

average_berX=sum(berX_all)/length(berX_all);
average_berA=sum(berA_all)/length(berA_all);
average_berB=sum(berB_all)/length(berB_all);
ber_all=[berX_all' berA_all' berB_all'];
%save('ber_all.mat', 'ber_all');
evmA=phy_cal_evm_by_pilot(measureEVMA,Num_data_symbols);
evmB=phy_cal_evm_by_pilot(measureEVMB,Num_data_symbols);
