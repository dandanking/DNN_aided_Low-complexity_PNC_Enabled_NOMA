%% PNC Packet Decoder
%  @author: lzyou@inc.cuhk.edu.hk
%  @date:   Nov 14, 2013
%
% Packet Format
%   A: TS1A      TS2A       DATA_A  CRC32
%   B:      TS1B      TS2B  DATA_B  CRC32

function [XBits, ABits, BBits, enXBits, enABits, enBBits,measureA,measureB] = PNC_Packet_Decoder(rxdata, Num_data_symbols, txdata_a, txdata_b)
    % Input:
    %   rxdata: rx samples (including noise)
    %   Num_data_symbols: number of data symbols
    %   txdata_a/b: tx samples of user a/b
    % Output:
    %   XBits: decoded source XOR bits
    %   ABits: decoded source A bits
    %   BBits: decoded source B bits
    %   en*Bits: decoded * bits (before interleaver and channel decoding)
    
    fft_length = 64;
    data_tones = 48;
    
    Num_samples_data = 80*(Num_data_symbols+4);
   
    rxdata = rxdata(1:Num_samples_data);
    TS1A = rxdata(1:80);
    TS1B = rxdata(81:160);
    rxTS2A = rxdata(161:240); txTS2A = txdata_a(161:240);
    rxTS2B = rxdata(241:320); txTS2B = txdata_b(241:320);
    RawData = rxdata(321:end);
        
    %----------to compensat frequency offset in time domain-----------%
    %[Delta_Theta_a] = phy_cfo_estimation (TS1A);
    %[Delta_Theta_b] = phy_cfo_estimation (TS1B);
    %Delta = (Delta_Theta_a + Delta_Theta_b)/2;
    %RawData = RawData .* exp(-j * [1:length(RawData)]' * Delta);

    %---------- to cut CP -----------%
    CutPos = 6;
    
    %---------- for channel estimation -----------%
    % for old data like benchmark_data/pnc_pkt.16S.32f.1.dat
    % H_a = My_Channel_Est (TS2A, Delta_Theta_a, CutPos);
    % H_b = My_Channel_Est (TS2B, Delta_Theta_b, CutPos);

    % for standard
    %H_a = phy_channel_estimation (rxTS2A, txTS2A, CutPos);
    %H_b = phy_channel_estimation (rxTS2B, txTS2B, CutPos);
    
    % for pveviously collected data
    % amp(pre)=16, amp(data)=8. Hence we should normailize by 8 so that for pre=2, we get data=1
    global MOD
    if strcmp(MOD,'BPSK')
      norm=8;          % for BPSK
    elseif strcmp(MOD,'QPSK')
      norm=8/sqrt(2);  % for QPSK
    end
    txTS2A = txTS2A(16+1:16+fft_length);
    txTS2B = txTS2B(16+1:16+fft_length);
    TS2_freq_a = fftshift(fft(txTS2A))/norm;
    TS2_freq_b = fftshift(fft(txTS2B))/norm;

    rxTS2A = rxTS2A(CutPos+1:CutPos+fft_length);
    rxTS2B = rxTS2B(CutPos+1:CutPos+fft_length);
    rxTS2_freq_a = fftshift(fft(rxTS2A));
    rxTS2_freq_b = fftshift(fft(rxTS2B));
    
    H_a = rxTS2_freq_a ./ TS2_freq_a;
    H_b = rxTS2_freq_b ./ TS2_freq_b;

    nremain = mod(length(RawData),80);
    if nremain ~= 0
        RawData = [RawData; zeros(80-nremain,1)];
    end
    
    fftRawData = zeros(64/80*length(RawData),1);
    m = 1;
    for i = 1 : 80 : length(RawData)
        timeRawData = RawData(i+CutPos:i+CutPos+fft_length-1);
        freqRawData = fft(timeRawData);
        freqRawData = fftshift(freqRawData);
        fftRawData((m-1)*fft_length+1 : m*fft_length) = freqRawData;
        m = m + 1;
    end

    null_encoded_bits=-10*ones(Num_data_symbols*data_tones,1);
    knownTxRawBits = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
    [RawBits, debug, measureA, measureB] = PNC_Data_Decoder ( fftRawData, data_tones, Num_data_symbols, H_a, H_b, knownTxRawBits);    
    

    if strcmp(MOD,'BPSK')
        A_INDEX=1; B_INDEX=2; X_INDEX=3;
        
        enXBits = RawBits(:,X_INDEX);
        enABits = RawBits(:,A_INDEX);
        enBBits = RawBits(:,B_INDEX);
    elseif strcmp(MOD,'QPSK')
        A_INDEX=[1 2]; B_INDEX=[3 4]; X_INDEX=[5 6];
        ABits = RawBits(:,A_INDEX); BBits = RawBits(:,B_INDEX); XBits = RawBits(:,X_INDEX); 
        enABits = reshape(ABits.', length(ABits)*2, 1); 
        enBBits = reshape(BBits.', length(BBits)*2, 1);
        enXBits = reshape(XBits.', length(XBits)*2, 1);
    end
        
    XBitsDeit = phy_de_interleaver(enXBits);
    ABitsDeit = phy_de_interleaver(enABits);
    BBitsDeit = phy_de_interleaver(enBBits);
    
    global DECODER_TYPE
    XBits = phy_viterbi_decoder(XBitsDeit,DECODER_TYPE);
    ABits = phy_viterbi_decoder(ABitsDeit,DECODER_TYPE); 
    BBits = phy_viterbi_decoder(BBitsDeit,DECODER_TYPE);
    
    %remove padding bits
    npadding = 8*1;
    XBits = XBits(1:end-npadding);
    ABits = ABits(1:end-npadding);
    BBits = BBits(1:end-npadding);
    

    
    
    
    