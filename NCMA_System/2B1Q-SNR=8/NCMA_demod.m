%  author: ph014@ie.cuhk.edu.hk
%  date:   April 20, 2015 / Modified 2020/11/16
%  demodulatoin

function [CiBits, CqBits, ABits, BBits, ABBits, ACiBits, ACqBits, BCiBits, BCqBits, ABCiBits, ABCqBits] ...
    = NCMA_demod(rx, hci, hcq, ha, hb, modulation, pkt_length, snr )

logp_A = zeros(pkt_length,1);
p_diff_A = zeros(pkt_length,1);
logp_B = zeros(pkt_length,1);
p_diff_B = zeros(pkt_length,1);
logp_Ci = zeros(pkt_length,1);
p_diff_Ci = zeros(pkt_length,1);
logp_Cq = zeros(pkt_length,1);
p_diff_Cq = zeros(pkt_length,1);

logp_xorAB = zeros(pkt_length,1);
p_diff_xorAB = zeros(pkt_length,1);
logp_xorACi = zeros(pkt_length,1);
p_diff_xorACi = zeros(pkt_length,1);
logp_xorACq = zeros(pkt_length,1);
p_diff_xorACq = zeros(pkt_length,1);
logp_xorBCi = zeros(pkt_length,1);
p_diff_xorBCi = zeros(pkt_length,1);
logp_xorBCq = zeros(pkt_length,1);
p_diff_xorBCq = zeros(pkt_length,1);
logp_xorABCi = zeros(pkt_length,1);
p_diff_xorABCi = zeros(pkt_length,1);
logp_xorABCq = zeros(pkt_length,1);
p_diff_xorABCq = zeros(pkt_length,1);

pA_0 = zeros(pkt_length,1);
pA_1 = zeros(pkt_length,1);
pB_0 = zeros(pkt_length,1);
pB_1 = zeros(pkt_length,1);
pCi_0 = zeros(pkt_length,1);
pCi_1 = zeros(pkt_length,1);
pCq_0 = zeros(pkt_length,1);
pCq_1 = zeros(pkt_length,1);
pXAB_0 = zeros(pkt_length,1);
pXAB_1 = zeros(pkt_length,1);
pXACi_0 = zeros(pkt_length,1);
pXACi_1 = zeros(pkt_length,1);
pXACq_0 = zeros(pkt_length,1);
pXACq_1 = zeros(pkt_length,1);
pXBCi_0 = zeros(pkt_length,1);
pXBCi_1 = zeros(pkt_length,1);
pXBCq_0 = zeros(pkt_length,1);
pXBCq_1 = zeros(pkt_length,1);
pXABCi_0 = zeros(pkt_length,1);
pXABCi_1 = zeros(pkt_length,1);
pXABCq_0 = zeros(pkt_length,1);
pXABCq_1 = zeros(pkt_length,1);

%compute probability
for ii=1:length(rx)
    snr_lin=10^(-snr/10);
    if strcmp(modulation,'BPSK') == 1
        data_pair = [1 1 1 1;
            1 -1 1 1;
            -1 1 1 1;
            -1 -1 1 1;
            1 1 -1 1;
            1 -1 -1 1;
            -1 1 -1 1;
            -1 -1 -1 1;
            1 1 1 -1;
            1 -1 1 -1;
            -1 1 1 -1;
            -1 -1 1 -1;
            1 1 -1 -1;
            1 -1 -1 -1;
            -1 1 -1 -1;
            -1 -1 -1 -1;]; 
    end
    
    d = zeros(length(data_pair),1);
    for pair_index = 1 : length(data_pair)
        d(pair_index) = norm( hci(ii) * data_pair(pair_index, 1)  ... 
            + hcq(ii) * data_pair(pair_index,2)  ...
            + ha(ii) * data_pair(pair_index,3)  ...
            + hb(ii) * data_pair(pair_index,4)  ...
            - rx(ii) )^2;
    end
    
    if strcmp(modulation,'BPSK') == 1
        
        %packet Ci
        p0_Ci = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_Ci = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_Ci(ii) = p0_Ci-p1_Ci;
        logp_Ci(ii) = log(p0_Ci-p1_Ci);
        pCi_0(ii) = p0_Ci;
        pCi_1(ii) = p1_Ci;
        
        %packet Cq
        p0_Cq = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin)  ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_Cq = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_Cq(ii) = p0_Cq-p1_Cq;
        logp_Cq(ii) = log(p0_Cq-p1_Cq);
        pCq_0(ii) = p0_Cq;
        pCq_1(ii) = p1_Cq;
        
        %packet A
        p0_A = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_A = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin);
        p_diff_A(ii) = p0_A-p1_A;
        logp_A(ii) = log(p0_A-p1_A);
        pA_0(ii) = p0_A;
        pA_1(ii) = p1_A;
        
        %packet B
        p0_B = exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_B = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p_diff_B(ii) = p0_B-p1_B;
        logp_B(ii) = log(p0_B-p1_B);
        pB_0(ii) = p0_B;
        pB_1(ii) = p1_B;
        
        %packet AxorB
        p0_xorAB = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorAB = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin);
        p_diff_xorAB(ii) = p0_xorAB-p1_xorAB;
        logp_xorAB(ii) = log(p0_xorAB-p1_xorAB);
        pXAB_0(ii) = p0_xorAB;
        pXAB_1(ii) = p1_xorAB;
        
        %packet AxorCi
        p0_xorACi = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorACi = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xorACi(ii) = p0_xorACi-p1_xorACi;
        logp_xorACi(ii) = log(p0_xorACi-p1_xorACi);
        pXACi_0(ii) = p0_xorACi;
        pXACi_1(ii) = p1_xorACi;
        
        %packet AxorCq
        p0_xorACq = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorACq = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xorACq(ii) = p0_xorACq-p1_xorACq;
        logp_xorACq(ii) = log(p0_xorACq-p1_xorACq);
        pXACq_0(ii) = p0_xorACq;
        pXACq_1(ii) = p1_xorACq;
        
        %packet BxorCi
        p0_xorBCi = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorBCi = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xorBCi(ii) = p0_xorBCi-p1_xorBCi;
        logp_xorBCi(ii) = log(p0_xorBCi-p1_xorBCi);
        pXBCi_0(ii) = p0_xorBCi;
        pXBCi_1(ii) = p1_xorBCi;
        
        %packet BxorCq
        p0_xorBCq = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorBCq = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xorBCq(ii) = p0_xorBCq-p1_xorBCq;
        logp_xorBCq(ii) = log(p0_xorBCq-p1_xorBCq);
        pXBCq_0(ii) = p0_xorBCq;
        pXBCq_1(ii) = p1_xorBCq;
        
        %packet AxorBxorCi
        p0_xorABCi = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorABCi = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xorABCi(ii) = p0_xorABCi-p1_xorABCi;
        logp_xorABCi(ii) = log(p0_xorABCi-p1_xorABCi);
        pXABCi_0(ii) = p0_xorABCi;
        pXABCi_1(ii) = p1_xorABCi;
        
        %packet AxorBxorCq
        p0_xorABCq = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xorABCq = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xorABCq(ii) = p0_xorABCq-p1_xorABCq;
        logp_xorABCq(ii) = log(p0_xorABCq-p1_xorABCq);
        pXABCq_0(ii) = p0_xorABCq;
        pXABCq_1(ii) = p1_xorABCq;
    end
    CiBits = -1*reshape([pCi_0,pCi_1]',1,length(pCi_0)*2);
    CqBits = -1*reshape([pCq_0,pCq_1]',1,length(pCq_0)*2);
    ABits = -1*reshape([pA_0,pA_1]',1,length(pA_0)*2); 
    BBits = -1*reshape([pB_0,pB_1]',1,length(pB_0)*2);
    ABBits = -1*reshape([pXAB_0,pXAB_1]',1,length(pXAB_0)*2);
    ACiBits = -1*reshape([pXACi_0,pXACi_1]',1,length(pXACi_0)*2);
    ACqBits = -1*reshape([pXACq_0,pXACq_1]',1,length(pXACq_0)*2);
    BCiBits = -1*reshape([pXBCi_0,pXBCi_1]',1,length(pXBCi_0)*2);
    BCqBits = -1*reshape([pXBCq_0,pXBCq_1]',1,length(pXBCq_0)*2);
    ABCiBits = -1*reshape([pXABCi_0,pXABCi_1]',1,length(pXABCi_0)*2);
    ABCqBits = -1*reshape([pXABCq_0,pXABCq_1]',1,length(pXABCq_0)*2);
end
end

