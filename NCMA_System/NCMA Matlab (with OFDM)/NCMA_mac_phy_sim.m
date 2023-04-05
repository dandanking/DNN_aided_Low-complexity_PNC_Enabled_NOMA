%% NCMA Trace-Driven Simulations: Cross-Layer Part
%  (1) PHY simulation (including two-phase decoding) is implemented elsewhere;
%  (2) Here we consider MAC layer throughput, and cross-layer simulations;
%
%% Mode:
%  (1) MUD  : MUD
%  (2) NCMA-: MUD + PHY bridging
%  (3) NCMA : MUD + PHY/MAC bridging
%  (4) NCMA+: MUD + PHY/MAC bridging + PHY SIC
%  (5) NCMA#: MUD + PHY/MAC bridging + PHY SIC + Cross-Layer forward decoding
%  (6) NCMA@: MUD + PHY/MAC bridging + PHY SIC + Cross-layer backward decoding
%  (7) NCMA*: MUD + PHY/MAC bridging + PHY SIC + Cross-layer forward/backward decoding
%
%% Input:
%  (1) phyRawMap: MUD map
%  (2) phyBridgeMap: PHY bridging map
%  (3) phyTwoPhaseMap: PHY bridging map + two-phase decoding map
%
%  author: lzyou@inc.cuhk.edu.hk
%  date:   April 20, 2013

clear;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
%-tx data
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
%tx_encoded_bits=[txA_encoded_bits,txB_encoded_bits,txX_encoded_bits];

txA = read_complex_binary('ncma_txdata1/txA_64S.dat');
txA_TS_time = txA(161+16:160+80);
TS2_freq_a = fftshift(fft(txA_TS_time))/8;    
txB = read_complex_binary('ncma_txdata1/txB_64S.dat');
txB_TS_time = txB(241+16:240+80);
TS2_freq_b = fftshift(fft(txB_TS_time))/8;


global DECODER_TYPE
DECODER_TYPE='soft';

fft_size=64;
data_tones=48;
padding_bytes=1;
num_data_sym=64;
nbits=num_data_sym*data_tones;

NULL_INDEX = [1:6 33 60:64];
A_INDEX = 1;
B_INDEX = 2;
X_INDEX = 3;

XOR_RECOVER_SIGNAL = -10;

% MODE: indicating which mode
NCMA_MODE = 1; SIC_MODE = 2; FULL_MODE = 3;
MODE([NCMA_MODE,SIC_MODE,FULL_MODE]) = [0,0,1];

% ------------------------------------------------------------------------%
%                          PHY Parameters                                 %
% ------------------------------------------------------------------------%

ncma_mac_files={'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
%ncma_mac_files={'0126_09_8dB'};
%ncma_mac_files={'0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
%ncma_mac_files={'0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB'};
%ncma_mac_files={'0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB'};
%ncma_mac_files={'0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB','0131_10_74dB_123dB'};
%ncma_mac_files={'0131_08_01_70dB_90dB'};
%ncma_mac_files={'75dB_75dB_rayleigh','75dB_95dB_rayleigh','75dB_105dB_rayleigh','75dB_135dB_rayleigh'};
%ncma_mac_files={'95dB_75dB_rayleigh','95dB_85dB_rayleigh','95dB_95dB_rayleigh','95dB_105dB_rayleigh','95dB_135dB_rayleigh'};
%ncma_mac_files={'5dB_awgn_random','6dB_awgn_random','7dB_awgn_random','8dB_awgn_random','9dB_awgn_random'};
%ncma_mac_files={'75dB_75dB_awgn_random','75dB_95dB_awgn_random','75dB_105dB_awgn_random','75dB_135dB_awgn_random'};
%ncma_mac_files={'95dB_75dB_awgn_random','95dB_85dB_awgn_random','95dB_95dB_awgn_random','95dB_105dB_awgn_random','95dB_135dB_awgn_random'};

for fileIndex=1:length(ncma_mac_files)

ncma_mac_tag=ncma_mac_files{fileIndex};

if MODE(NCMA_MODE)
    ss=strcat('matlab_data/rmud_phy/',ncma_mac_tag);
    ss=strcat(ss,'_maps.mat');
elseif MODE(SIC_MODE)
    ss=strcat('matlab_data/sic_phy/',ncma_mac_tag);
    ss =strcat(ss,'_softSIC_map.mat');
elseif MODE(FULL_MODE)
    ss=strcat('matlab_data/rmud_sic_phy/',ncma_mac_tag);
    ss=strcat(ss,'_softFull_map.mat');
end
load(ss);
if MODE(SIC_MODE)
    phyRawMap = softPhySIC2XorMap;
    phyBridgeMap = phy_gen_bridge_map(phyRawMap);
elseif MODE(FULL_MODE)
    phyBridgeMap = softPhyFullMap;
end

fprintf(' Opening file: %s \n', ss);

%--------------------------------------------%
ss=strcat('ncma_demod_data/',ncma_mac_tag);
%ss=strcat('ncma_demod_data/rayleigh_unbalanced_75dB/',ncma_mac_tag);
%ss=strcat('ncma_demod_data/rayleigh_unbalanced_95dB/',ncma_mac_tag);
%ss=strcat('ncma_demod_data/awgn_random_degree/',ncma_mac_tag);
%ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_75dB/',ncma_mac_tag);
%ss=strcat('ncma_demod_data/awgn_random_degree_unbalanced_95dB/',ncma_mac_tag);
%--------------------------------------------%
ss=strcat(ss,'/rx-ofdmdemod.dat');
demod_in = read_complex_binary(ss);

ss=strcat(ss,'b');
demod_flag = read_char_binary(ss);

global fphy
global fmac
ss1=strcat(ncma_mac_tag,'_phy.log'); ss2=strcat(ncma_mac_tag,'_mac.log');
%fphy = fopen(ss1,'wb'); fmac = fopen(ss2,'wb');
fphy = 1; fmac = 1;

pos_list=find(demod_flag>0);
num_pkts=length(pos_list);

% ------------------------------------------------------------------------%
%                          MAC Parameters                                 %
% ------------------------------------------------------------------------%

RANGE_L = [4,8,16,32];
RANGE_L = [16];
nL = length(RANGE_L);

RANGE_RATIO = [0.5,0.666,1,1.5,2];
%RANGE_RATIO = [1,1.5,2,3,4,5];
RANGE_RATIO = [2,3];
nratio = length(RANGE_RATIO);


% Format: MUD, NCMA-, NCMA, NCMA+, NCMA#, NCMA@, NCMA*
INDEX={1,2,3,4,5,6,7,8,9,10,11,12,13,14};
[A_MUD_INDEX,B_MUD_INDEX,A_NCMA1_INDEX,B_NCMA1_INDEX, ...
 A_NCMA2_INDEX,B_NCMA2_INDEX,A_NCMA3_INDEX,B_NCMA3_INDEX, ...
 A_NCMA4_INDEX,B_NCMA4_INDEX,A_NCMA5_INDEX,B_NCMA5_INDEX, ...
 A_NCMA6_INDEX,B_NCMA6_INDEX] = deal(INDEX{:});

npktsAll = zeros(nratio*nL,length(INDEX));

% Format: MUD, NCMA-, NCMA, NCMA_UPPER, NCMA+, NCMA+_UPPER, NCMA#, NCMA#_UPPER, NCMA@, NCMA@_UPPER, NCMA*, NCMA*_UPPER
INDEX={1,2,3,4,5,6,7,8,9,10,11,12};
[MUD_INDEX, NCMA_MINUS_INDEX, NCMA_INDEX, NCMA_UPPER_INDEX, ...
 NCMA_PLUS_INDEX, NCMA_PLUS_UPPER_INDEX, NCMA_SHOP_INDEX, NCMA_SHOP_UPPER_INDEX, ...
 NCMA_AT_INDEX, NCMA_AT_UPPER_INDEX, NCMA_STAR_INDEX, NCMA_STAR_UPPER_INDEX] = deal(INDEX{:});

npkts = zeros(nratio*nL,length(INDEX));
SLOTS = 1000;

% FLAG: indicating which mode
FLAG = zeros(length(INDEX),1);
if MODE(NCMA_MODE)
    FLAG([MUD_INDEX, NCMA_MINUS_INDEX, NCMA_INDEX, NCMA_PLUS_INDEX, NCMA_SHOP_INDEX, NCMA_AT_INDEX, NCMA_STAR_INDEX]) = [1,0,0,0,0,0,0]; %[1,1,1,1,1,1,1];
    postfix = '_mac.mat';
elseif MODE(SIC_MODE)
    FLAG([MUD_INDEX, NCMA_MINUS_INDEX, NCMA_INDEX, NCMA_PLUS_INDEX, NCMA_SHOP_INDEX, NCMA_AT_INDEX, NCMA_STAR_INDEX]) = [1,1,1,0,0,0,0];
    postfix = '_sic_mac.mat';
elseif MODE(FULL_MODE)
    FLAG([MUD_INDEX, NCMA_MINUS_INDEX, NCMA_INDEX, NCMA_PLUS_INDEX, NCMA_SHOP_INDEX, NCMA_AT_INDEX, NCMA_STAR_INDEX]) = [0,0,1,0,0,0,0];
    postfix = '_full_mac.mat';
end
% ------------------------------------------------------------------------%
%                           Main Functions                                %
% ------------------------------------------------------------------------%
for rr=1:length(RANGE_RATIO)
 
ratio = RANGE_RATIO(rr);

for ll=1:length(RANGE_L)
    
L = RANGE_L(ll);

L_A = round(L*ratio);
L_B = L;
L_X = max(L_A, L_B);

fprintf(fmac,' rr=%d ll=%d L=%d ratio=%f\n', rr, ll, L, ratio);


d_ncma2_X_decodable = 0;
d_ncma3_X_decodable = 0;
d_ncma4_X_decodable = 0;
d_ncma5_X_decodable = 0;
d_ncma6_X_decodable = 0;

d_ncma3_X_cnt = 0;

%FIXME: implement detection map, and single-user PHY decoder
%phyDetectMap = zeros(num_pkts,3);       % packet detection map

mudMACMap   = zeros(num_pkts,3);        % MUD  : mud decoding map
ncma1MACMap = zeros(num_pkts,3);        % NCMA-: mud decoding map + PHY bridging
ncma2MACMap = zeros(num_pkts,3);        % NCMA : mud decoding map + PHY/MAC bridging
ncma3MACMap = zeros(num_pkts,3);        % NCMA+: mud decoding map + PHY/MAC bridging + PHY SIC
ncma4MACMap = zeros(num_pkts,3);        % NCMA#: mud decoding map + PHY/MAC bridging + PHY SIC + Cross-Layer forward decoding
ncma5MACMap = zeros(num_pkts,3);        % NCMA@: mud decoding map + PHY/MAC bridging + PHY SIC + Cross-layer backward decoding
ncma6MACMap = zeros(num_pkts,3);        % NCMA*: mud decoding map + PHY/MAC bridging + PHY/Cross-layer forward/backward decoding

% ncma1: NCMA-; ncma2: NCMA; ncma3: NCMA+; ncma4: NCMA#; ncma5: NCMA@; ncma6: NCMA*
mudCntA = 0; mudCntB = 0;
ncma1CntA = 0; ncma1CntB = 0;
ncma2CntA = 0; ncma2CntB = 0; ncma2CntX = 0;
ncma3CntA = 0; ncma3CntB = 0; ncma3CntX = 0;
ncma4CntA = 0; ncma4CntB = 0; ncma4CntX = 0;
ncma5CntA = 0; ncma5CntB = 0; ncma5CntX = 0;
ncma6CntA = 0; ncma6CntB = 0; ncma6CntX = 0;

mudPktNumA = 1; mudPktNumB = 1;
ncma1PktNumA = 1; ncma1PktNumB = 1;
ncma2PktNumA = 1; ncma2PktNumB = 1; ncma2PktNumX = 1;
ncma3PktNumA = 1; ncma3PktNumB = 1; ncma3PktNumX = 1;
ncma4PktNumA = 1; ncma4PktNumB = 1; ncma4PktNumX = 1;
ncma5PktNumA = 1; ncma5PktNumB = 1; ncma5PktNumX = 1;
ncma6PktNumA = 1; ncma6PktNumB = 1; ncma6PktNumX = 1;

% ------------------------------------------------------------------------%
%                            Main Function                                %
% ------------------------------------------------------------------------%
for pktIndex=1:num_pkts
    nsym = num_data_sym;
            
    pos = pos_list(pktIndex);
    rxdata = demod_in((pos-1)*fft_size+1:(pos+nsym-1+2)*fft_size);
    rxTS2A = rxdata(1:64);
    rxTS2B = rxdata(65:128);
    H_a = rxTS2A ./ TS2_freq_a;
    H_b = rxTS2B ./ TS2_freq_b;
    rxRawData = rxdata(129:end);

    %---------------------------------------------------------------------%
    %                         MUD MAC Layer                               %
    %---------------------------------------------------------------------%
    if FLAG(MUD_INDEX)
        okA=phyRawMap(pktIndex,A_INDEX);
        okB=phyRawMap(pktIndex,B_INDEX);
        if okA
            mudCntA = mudCntA + 1;
            mudMACMap(pktIndex,A_INDEX) = mudPktNumA;
        end
        if okB
            mudCntB = mudCntB + 1;
            mudMACMap(pktIndex,B_INDEX) = mudPktNumB;
        end
        
        if mudCntA == L_A
            mudPktNumA = mudPktNumA + 1;
            mudCntA = 0;
        end
        if mudCntB == L_B
            mudPktNumB = mudPktNumB + 1;
            mudCntB = 0;
        end
    end
    %---------------------------------------------------------------------%
    %                     NCMA-/NCMA MAC Layer                            %
    %---------------------------------------------------------------------%
    % NCMA-(ncma1): only PHY bridging
    % NCMA (ncma2):  PHY/MAC bridging
    %---------------------------------------------------------------------%
    if FLAG(NCMA_MINUS_INDEX) || FLAG(NCMA_INDEX)
        okA=phyBridgeMap(pktIndex,A_INDEX);
        okB=phyBridgeMap(pktIndex,B_INDEX);
        okX=phyBridgeMap(pktIndex,X_INDEX);
        
        if okA && okB && ~okX
            okX = 1;
        elseif ~okA && okB && okX
            okA = 1;
        elseif okA && ~okB && okX
            okB = 1;
        end
    
        if FLAG(NCMA_MINUS_INDEX)
            if okA
                ncma1CntA = ncma1CntA + 1;
                ncma1MACMap(pktIndex,A_INDEX) = ncma1PktNumA;
            end
            if okB
                ncma1CntB = ncma1CntB + 1;
                ncma1MACMap(pktIndex,B_INDEX) = ncma1PktNumB;
            end
            
            if ncma1CntA == L_A
                ncma1PktNumA = ncma1PktNumA + 1;
                ncma1CntA = 0;
            end
            if ncma1CntB == L_B
                ncma1PktNumB = ncma1PktNumB + 1;
                ncma1CntB = 0;
            end
        end
        
        if FLAG(NCMA_INDEX)
            if d_ncma2_X_decodable && ~okX
                okX = 1;
                if ~okA && okB
                    okA = 1;
                elseif okA && ~okB
                    okB = 1;
                end
            end
            
            fprintf(fmac,'ncma2: curIndex=%d ok=[%d %d %d] countA=[%d %d] countB=[%d %d] countX=[%d %d] \n',pktIndex,okA,okB,okX,ncma2PktNumA,ncma2CntA,ncma2PktNumB,ncma2CntB,ncma2PktNumX,ncma2CntX);
            if okA
                ncma2CntA = ncma2CntA + 1;
                ncma2MACMap(pktIndex,A_INDEX) = ncma2PktNumA;
            end
            if okB
                ncma2CntB = ncma2CntB + 1;
                ncma2MACMap(pktIndex,B_INDEX) = ncma2PktNumB;
            end
            if okX
                ncma2CntX = ncma2CntX + 1;
                ncma2MACMap(pktIndex,X_INDEX) = ncma2PktNumX;
            end
            
            % NCMA: complicate than NCMA- since it has XOR equation and MAC bridging
            inputSeqInfo = zeros(3,2);
            inputSeqInfo(A_INDEX,:) = [ncma2CntA, ncma2PktNumA];
            inputSeqInfo(B_INDEX,:) = [ncma2CntB, ncma2PktNumB];
            inputSeqInfo(X_INDEX,:) = [ncma2CntX, ncma2PktNumX];
            %fprintf(' pktIndex=%d \n', pktIndex);
            if pktIndex == 77
                pktIndex = 77;
            end
            [ret,ncma2MACMap,d_ncma2_X_decodable] = mac_update_seqno(pktIndex,inputSeqInfo,d_ncma2_X_decodable,ncma2MACMap,L_A,L_B,'ncma2',...
                [],[],[],[],[],[],[],[],[],[],[]);
            ncma2CntA = ret(A_INDEX,1); ncma2PktNumA = ret(A_INDEX,2);
            ncma2CntB = ret(B_INDEX,1); ncma2PktNumB = ret(B_INDEX,2);
            ncma2CntX = ret(X_INDEX,1); ncma2PktNumX = ret(X_INDEX,2);
        end
    end
    %---------------------------------------------------------------------%
    %                     NCMA+/NCMA* MAC Layer                           %
    %---------------------------------------------------------------------%
    % NCMA+(ncma3): NCMA add PHY SIC
    % NCMA#(ncma4): add Cross-Layer forward decoding
    % NCMA@(ncma5): add Cross-Layer backward decoding
    % NCMA*(ncma6): add Cross-Layer forward/backward decoding
    %---------------------------------------------------------------------%   
    %FIXME: here we assume FLAG(NCMA+, NCMA#, NCMA*) are [1,1,1] or [0,0,0]
    if FLAG(NCMA_PLUS_INDEX) || FLAG(NCMA_SHOP_INDEX) || FLAG(NCMA_AT_INDEX) || FLAG(NCMA_STAR_INDEX)
        okA=phyTwoPhaseMap(pktIndex,A_INDEX);
        okB=phyTwoPhaseMap(pktIndex,B_INDEX);
        okX=phyTwoPhaseMap(pktIndex,X_INDEX);
        
        % In case we do not do Bridging
        if okA && okB && ~okX
            okX = 1;
        elseif ~okA && okB && okX
            okA = 1;
        elseif okA && ~okB && okX
            okB = 1;
        end
        
        % for NCMA+/NCMA#/NCMA@/NCMA*
        if okA
            if FLAG(NCMA_PLUS_INDEX)
                ncma3CntA = ncma3CntA + 1;
                ncma3MACMap(pktIndex,A_INDEX) = ncma3PktNumA;
            end
            if FLAG(NCMA_SHOP_INDEX)
                ncma4CntA = ncma4CntA + 1;
                ncma4MACMap(pktIndex,A_INDEX) = ncma4PktNumA;
            end
            if FLAG(NCMA_AT_INDEX)
                ncma5CntA = ncma5CntA + 1;
                ncma5MACMap(pktIndex,A_INDEX) = ncma5PktNumA;
            end
            if FLAG(NCMA_STAR_INDEX)
                ncma6CntA = ncma6CntA + 1;
                ncma6MACMap(pktIndex,A_INDEX) = ncma6PktNumA;
            end
        end
        if okB
            if FLAG(NCMA_PLUS_INDEX)
                ncma3CntB = ncma3CntB + 1;
                ncma3MACMap(pktIndex,B_INDEX) = ncma3PktNumB;
            end
            if FLAG(NCMA_SHOP_INDEX)
                ncma4CntB = ncma4CntB + 1;
                ncma4MACMap(pktIndex,B_INDEX) = ncma4PktNumB;
            end
            if FLAG(NCMA_AT_INDEX)
                ncma5CntB = ncma5CntB + 1;
                ncma5MACMap(pktIndex,B_INDEX) = ncma5PktNumB;
            end
            if FLAG(NCMA_STAR_INDEX)
                ncma6CntB = ncma6CntB + 1;
                ncma6MACMap(pktIndex,B_INDEX) = ncma6PktNumB;
            end
        end
        if okX
            if FLAG(NCMA_PLUS_INDEX)
                ncma3CntX = ncma3CntX + 1;
                ncma3MACMap(pktIndex,X_INDEX) = ncma3PktNumX;
            end
            if FLAG(NCMA_SHOP_INDEX)
                ncma4CntX = ncma4CntX + 1;
                ncma4MACMap(pktIndex,X_INDEX) = ncma4PktNumX;
            end
            if FLAG(NCMA_AT_INDEX)
                ncma5CntX = ncma5CntX + 1;
                ncma5MACMap(pktIndex,X_INDEX) = ncma5PktNumX;
            end
            if FLAG(NCMA_STAR_INDEX)
                ncma6CntX = ncma6CntX + 1;
                ncma6MACMap(pktIndex,X_INDEX) = ncma6PktNumX;
            end
        end
        
        % If we have recovered X, we can apply Cross-Layer forward decoding
        % No need to consider A/B here
        if ~okA && ~okB && ~okX
            knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
            if d_ncma4_X_decodable || d_ncma6_X_decodable
                berX = -1;
                knownTxRawData(:,X_INDEX) = txX_encoded_bits;
                RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                rawABits = RawBits(:,A_INDEX); rawBBits = RawBits(:,B_INDEX);
                ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),'soft');
                BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),'soft');
                ABits=ABits(1:end-8);
                BBits=BBits(1:end-8);
                [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits,txA_source_bits,0);
                [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits,txB_source_bits,0);
                if okA || okB
                    % only for NCMA4(NCMA#) and NCMA6(NCMA*)
                    if d_ncma4_X_decodable && FLAG(NCMA_SHOP_INDEX)
                        ncma4CntA = ncma4CntA + 1;
                        ncma4MACMap(pktIndex,A_INDEX) = ncma4PktNumA; %XOR_RECOVER_SIGNAL*ncma4PktNumA;
                        ncma4CntB = ncma4CntB + 1;
                        ncma4MACMap(pktIndex,B_INDEX) = ncma4PktNumB; %XOR_RECOVER_SIGNAL*ncma4PktNumB;
                    end
                    if d_ncma6_X_decodable && FLAG(NCMA_STAR_INDEX)
                        ncma6CntA = ncma6CntA + 1;
                        ncma6MACMap(pktIndex,A_INDEX) = ncma6PktNumA; %XOR_RECOVER_SIGNAL*ncma5PktNumA;
                        ncma6CntB = ncma6CntB + 1;
                        ncma6MACMap(pktIndex,B_INDEX) = ncma6PktNumB; %XOR_RECOVER_SIGNAL*ncma5PktNumB;
                    end
                end
                
                if d_ncma4_X_decodable && FLAG(NCMA_SHOP_INDEX)
                    fprintf(fphy,'ncma4:  forward X+%d: okA=%d berA=%f | okB=%d berB=%f | countA=[%d %d] countB=[%d %d] countX=[%d %d] \n', pktIndex, okA, berA, okB, berB, ncma4PktNumA, ncma4CntA, ncma4PktNumB, ncma4CntB, ncma4PktNumX, ncma4CntX);
                end
                if d_ncma6_X_decodable && FLAG(NCMA_STAR_INDEX)
                    fprintf(fphy,'ncma6:  forward X+%d: okA=%d berA=%f | okB=%d berB=%f | countA=[%d %d] countB=[%d %d] countX=[%d %d] \n', pktIndex, okA, berA, okB, berB, ncma6PktNumA, ncma6CntA, ncma6PktNumB, ncma6CntB, ncma6PktNumX, ncma6CntX);
                end
            end
        end
        
        if FLAG(NCMA_PLUS_INDEX)
            % NCMA+
            inputSeqInfo = zeros(3,2);
            inputSeqInfo(A_INDEX,:) = [ncma3CntA, ncma3PktNumA];
            inputSeqInfo(B_INDEX,:) = [ncma3CntB, ncma3PktNumB];
            inputSeqInfo(X_INDEX,:) = [ncma3CntX, ncma3PktNumX];
            %fprintf(' pktIndex=%d countX=%d decodable=%d \n',pktIndex,ncma3CntX,d_ncma3_X_decodable);
            [ret,ncma3MACMap,d_ncma3_X_decodable] = mac_update_seqno(pktIndex,inputSeqInfo,d_ncma3_X_decodable,ncma3MACMap,L_A,L_B,'ncma3',...
                [],[],[],[],[],[],[],[],[],[],[]);
            ncma3CntA = ret(A_INDEX,1); ncma3PktNumA = ret(A_INDEX,2);
            ncma3CntB = ret(B_INDEX,1); ncma3PktNumB = ret(B_INDEX,2);
            ncma3CntX = ret(X_INDEX,1); ncma3PktNumX = ret(X_INDEX,2);
        end
        
        if FLAG(NCMA_SHOP_INDEX)
            % NCMA#
            inputSeqInfo = zeros(3,2);
            inputSeqInfo(A_INDEX,:) = [ncma4CntA, ncma4PktNumA];
            inputSeqInfo(B_INDEX,:) = [ncma4CntB, ncma4PktNumB];
            inputSeqInfo(X_INDEX,:) = [ncma4CntX, ncma4PktNumX];
            %fprintf(' pktIndex=%d countX=%d decodable=%d \n',pktIndex,ncma4CntX,d_ncma4_X_decodable);
            [ret,ncma4MACMap,d_ncma4_X_decodable] = mac_update_seqno(pktIndex,inputSeqInfo,d_ncma4_X_decodable,ncma4MACMap,L_A,L_B,'ncma4',...
                [],[],[],[],[],[],[],[],[],[],[]);
            ncma4CntA = ret(A_INDEX,1); ncma4PktNumA = ret(A_INDEX,2);
            ncma4CntB = ret(B_INDEX,1); ncma4PktNumB = ret(B_INDEX,2);
            ncma4CntX = ret(X_INDEX,1); ncma4PktNumX = ret(X_INDEX,2);
        end

        if FLAG(NCMA_AT_INDEX)
            % NCMA@ (backward decoding)
            inputSeqInfo = zeros(3,2);
            inputSeqInfo(A_INDEX,:) = [ncma5CntA, ncma5PktNumA];
            inputSeqInfo(B_INDEX,:) = [ncma5CntB, ncma5PktNumB];
            inputSeqInfo(X_INDEX,:) = [ncma5CntX, ncma5PktNumX];
            %fprintf('ncma5: pktIndex=%d countA=[%d %d] countB=[%d %d] countX=[%d %d] \n',pktIndex,ncma5PktNumA,ncma5CntA,ncma5PktNumB,ncma5CntB,ncma5PktNumX,ncma5CntX);
            [ret,ncma5MACMap,d_ncma5_X_decodable] = mac_update_seqno(pktIndex,inputSeqInfo,d_ncma5_X_decodable,ncma5MACMap,L_A,L_B, 'ncma5',...
                rxRawData,data_tones,nsym,H_a,H_b, ...
                txA_encoded_bits,txB_encoded_bits,txX_encoded_bits, ...
                txA_source_bits,txB_source_bits,txX_source_bits);
            ncma5CntA = ret(A_INDEX,1); ncma5PktNumA = ret(A_INDEX,2);
            ncma5CntB = ret(B_INDEX,1); ncma5PktNumB = ret(B_INDEX,2);
            ncma5CntX = ret(X_INDEX,1); ncma5PktNumX = ret(X_INDEX,2);
        end
        
        if FLAG(NCMA_STAR_INDEX)
            % NCMA* (backward decoding)
            inputSeqInfo = zeros(3,2);
            inputSeqInfo(A_INDEX,:) = [ncma6CntA, ncma6PktNumA];
            inputSeqInfo(B_INDEX,:) = [ncma6CntB, ncma6PktNumB];
            inputSeqInfo(X_INDEX,:) = [ncma6CntX, ncma6PktNumX];
            %fprintf('ncma6: pktIndex=%d countA=[%d %d] countB=[%d %d] countX=[%d %d] decodable=%d \n',pktIndex,ncma6PktNumA,ncma6CntA,ncma6PktNumB,ncma6CntB,ncma6PktNumX,ncma6CntX,d_ncma6_X_decodable);
            [ret,ncma6MACMap,d_ncma6_X_decodable] = mac_update_seqno(pktIndex,inputSeqInfo,d_ncma6_X_decodable,ncma6MACMap,L_A,L_B, 'ncma6',...
                rxRawData,data_tones,nsym,H_a,H_b, ...
                txA_encoded_bits,txB_encoded_bits,txX_encoded_bits, ...
                txA_source_bits,txB_source_bits,txX_source_bits);
            ncma6CntA = ret(A_INDEX,1); ncma6PktNumA = ret(A_INDEX,2);
            ncma6CntB = ret(B_INDEX,1); ncma6PktNumB = ret(B_INDEX,2);
            ncma6CntX = ret(X_INDEX,1); ncma6PktNumX = ret(X_INDEX,2);
        end
    end
    %-------------------------------%
end

cur = (rr-1)*nL+ll;

npktsAll(cur,A_MUD_INDEX) = ((mudPktNumA-1)*L_A+mudCntA)/SLOTS;
npktsAll(cur,B_MUD_INDEX) = ((mudPktNumB-1)*L_B+mudCntB)/SLOTS;
npktsAll(cur,A_NCMA1_INDEX) = ((ncma1PktNumA-1)*L_A+ncma1CntA)/SLOTS;
npktsAll(cur,B_NCMA1_INDEX) = ((ncma1PktNumB-1)*L_B+ncma1CntB)/SLOTS;
npktsAll(cur,A_NCMA2_INDEX) = ((ncma2PktNumA-1)*L_A+ncma2CntA)/SLOTS;
npktsAll(cur,B_NCMA2_INDEX) = ((ncma2PktNumB-1)*L_B+ncma2CntB)/SLOTS;
npktsAll(cur,A_NCMA3_INDEX) = ((ncma3PktNumA-1)*L_A+ncma3CntA)/SLOTS;
npktsAll(cur,B_NCMA3_INDEX) = ((ncma3PktNumB-1)*L_B+ncma3CntB)/SLOTS;
npktsAll(cur,A_NCMA4_INDEX) = ((ncma4PktNumA-1)*L_A+ncma4CntA)/SLOTS;
npktsAll(cur,B_NCMA4_INDEX) = ((ncma4PktNumB-1)*L_B+ncma4CntB)/SLOTS;
npktsAll(cur,A_NCMA5_INDEX) = ((ncma5PktNumA-1)*L_A+ncma5CntA)/SLOTS;
npktsAll(cur,B_NCMA5_INDEX) = ((ncma5PktNumB-1)*L_B+ncma5CntB)/SLOTS;
npktsAll(cur,A_NCMA6_INDEX) = ((ncma6PktNumA-1)*L_A+ncma6CntA)/SLOTS;
npktsAll(cur,B_NCMA6_INDEX) = ((ncma6PktNumB-1)*L_B+ncma6CntB)/SLOTS;

npkts(cur,MUD_INDEX) = npktsAll(cur,A_MUD_INDEX)+npktsAll(cur,B_MUD_INDEX);
npkts(cur,NCMA_MINUS_INDEX) = npktsAll(cur,A_NCMA1_INDEX)+npktsAll(cur,B_NCMA1_INDEX);
npkts(cur,NCMA_INDEX) = npktsAll(cur,A_NCMA2_INDEX)+npktsAll(cur,B_NCMA2_INDEX);
npkts(cur,NCMA_PLUS_INDEX) = npktsAll(cur,A_NCMA3_INDEX)+npktsAll(cur,B_NCMA3_INDEX);
npkts(cur,NCMA_SHOP_INDEX) = npktsAll(cur,A_NCMA4_INDEX)+npktsAll(cur,B_NCMA4_INDEX);
npkts(cur,NCMA_AT_INDEX) = npktsAll(cur,A_NCMA5_INDEX)+npktsAll(cur,B_NCMA5_INDEX);
npkts(cur,NCMA_STAR_INDEX) = npktsAll(cur,A_NCMA6_INDEX)+npktsAll(cur,B_NCMA6_INDEX);

if FLAG(MUD_INDEX)
[two_equations_num,one_equations_num]=mac_count_equations(phyRawMap);
npkts(cur,NCMA_UPPER_INDEX) = (two_equations_num*2+one_equations_num)/SLOTS;
end

if FLAG(NCMA_PLUS_INDEX)
[two_equations_num,one_equations_num]=mac_count_equations(phyTwoPhaseMap);
npkts(cur,NCMA_PLUS_UPPER_INDEX) = (two_equations_num*2+one_equations_num)/SLOTS;
end

if FLAG(NCMA_SHOP_INDEX)
[two_equations_num,one_equations_num]=mac_count_equations(ncma4MACMap);
npkts(cur,NCMA_SHOP_UPPER_INDEX) = (two_equations_num*2+one_equations_num)/SLOTS;
end

if FLAG(NCMA_AT_INDEX)
[two_equations_num,one_equations_num]=mac_count_equations(ncma5MACMap);
npkts(cur,NCMA_AT_UPPER_INDEX) = (two_equations_num*2+one_equations_num)/SLOTS;
end

if FLAG(NCMA_STAR_INDEX)
[two_equations_num,one_equations_num]=mac_count_equations(ncma6MACMap);
npkts(cur,NCMA_STAR_UPPER_INDEX) = (two_equations_num*2+one_equations_num)/SLOTS;
end


npkts

end % end of RANGE_L

end % end of RANGE_RATIO

if fmac ~= 1
    fclose(fmac); 
end
if fphy ~= 1
    fclose(fphy);
end

clear demod_in;
clear demod_flag;

if MODE(NCMA_MODE) || MODE(FULL_MODE)
    ss=strcat('matlab_data/',ncma_mac_tag);
elseif MODE(SIC_MODE)
    ss=strcat('matlab_data/',ncma_mac_tag);
end
ss=strcat(ss,postfix);
save(ss);

end % end of FILES