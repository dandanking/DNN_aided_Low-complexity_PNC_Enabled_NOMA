%  author: ph014@ie.cuhk.edu.hk
%  date:   April 20, 2015 / Modified 2020/11/16
%  demodulatoin

function [ ABits, BBits, CBits, X1Bits, X2Bits, X3Bits, X0Bits, p_diff_A,p_diff_B,p_diff_xor1 ] ...
    = NCMA_demod(rx, ha, hb, hc, modulation, pkt_length, snr )

logp_A = zeros(pkt_length,1);
p_diff_A = zeros(pkt_length,1);
logp_B = zeros(pkt_length,1);
p_diff_B = zeros(pkt_length,1);
logp_C = zeros(pkt_length,1);
p_diff_C = zeros(pkt_length,1);

logp_xor1 = zeros(pkt_length,1);
p_diff_xor1 = zeros(pkt_length,1);
logp_xor2 = zeros(pkt_length,1);
p_diff_xor2 = zeros(pkt_length,1);
logp_xor3 = zeros(pkt_length,1);
p_diff_xor3 = zeros(pkt_length,1);
logp_xor0 = zeros(pkt_length,1);
p_diff_xor0 = zeros(pkt_length,1);

pA_0 = zeros(pkt_length,1);
pA_1 = zeros(pkt_length,1);
pB_0 = zeros(pkt_length,1);
pB_1 = zeros(pkt_length,1);
pC_0 = zeros(pkt_length,1);
pC_1 = zeros(pkt_length,1);
pX1_0 = zeros(pkt_length,1);
pX1_1 = zeros(pkt_length,1);
pX2_0 = zeros(pkt_length,1);
pX2_1 = zeros(pkt_length,1);
pX3_0 = zeros(pkt_length,1);
pX3_1 = zeros(pkt_length,1);
pX0_0 = zeros(pkt_length,1);
pX0_1 = zeros(pkt_length,1);

%compute probability
for ii=1:length(rx)
    snr_lin=10^(-snr/10);
    if strcmp(modulation,'BPSK') == 1
        data_pair = [1 1 1;
            1 -1 1;
            -1 1 1;
            -1 -1 1;
            1 1 -1;
            1 -1 -1;
            -1 1 -1;
            -1 -1 -1;];  
    end
    
    d = zeros(length(data_pair),1);
    for pair_index = 1 : length(data_pair)
        d(pair_index) = norm( ha(ii) * data_pair(pair_index, 1)  ... 
            + hb(ii) * data_pair(pair_index,2)  ...
            + hc(ii) * data_pair(pair_index,3)  ...
            - rx(ii) )^2;
    end
    
    if strcmp(modulation,'BPSK') == 1
        %packet A
        p0_A = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_A = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin);
        p_diff_A(ii) = p0_A-p1_A;
        logp_A(ii) = log(p0_A-p1_A);
        pA_0(ii) = p0_A;
        pA_1(ii) = p1_A;
        
        %packet B
        p0_B = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_B = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin);
        p_diff_B(ii) = p0_B-p1_B;
        logp_B(ii) = log(p0_B-p1_B);
        pB_0(ii) = p0_B;
        pB_1(ii) = p1_B;
        
        %packet C
        p0_C = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_C = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin);
        p_diff_C(ii) = p0_C-p1_C;
        logp_C(ii) = log(p0_C-p1_C);
        pC_0(ii) = p0_C;
        pC_1(ii) = p1_C;
        
        %packet xor1: AxorB
        p0_xor1 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_xor1 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin);
        p_diff_xor1(ii) = p0_xor1-p1_xor1;
        logp_xor1(ii) = log(p0_xor1-p1_xor1);
        pX1_0(ii) = p0_xor1;
        pX1_1(ii) = p1_xor1;
        
        %packet xor2: AxorC
        p0_xor2 = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_xor2 = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin);
        p_diff_xor2(ii) = p0_xor2-p1_xor2;
        logp_xor2(ii) = log(p0_xor2-p1_xor2);
        pX2_0(ii) = p0_xor2;
        pX2_1(ii) = p1_xor2;
        
        %packet xor3: BxorC
        p0_xor3 = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_xor3 = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin);
        p_diff_xor3(ii) = p0_xor3-p1_xor3;
        logp_xor3(ii) = log(p0_xor3-p1_xor3);
        pX3_0(ii) = p0_xor3;
        pX3_1(ii) = p1_xor3;
        
        %packet xor0: AxorBxorC
        p0_xor0 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p1_xor0 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin);
        p_diff_xor0(ii) = p0_xor0-p1_xor0;
        logp_xor0(ii) = log(p0_xor0-p1_xor0);
        pX0_0(ii) = p0_xor0;
        pX0_1(ii) = p1_xor0;
    end
    
    ABits = -1*reshape([pA_0,pA_1]',1,length(pA_0)*2); 
    %fprintf('packet = %d\n',ABits);
    BBits = -1*reshape([pB_0,pB_1]',1,length(pB_0)*2);
    CBits = -1*reshape([pC_0,pC_1]',1,length(pC_0)*2);
    X1Bits = -1*reshape([pX1_0,pX1_1]',1,length(pX1_0)*2);
    X2Bits = -1*reshape([pX2_0,pX2_1]',1,length(pX2_0)*2);
    X3Bits = -1*reshape([pX3_0,pX3_1]',1,length(pX3_0)*2);
    X0Bits = -1*reshape([pX0_0,pX0_1]',1,length(pX0_0)*2);
end
end

