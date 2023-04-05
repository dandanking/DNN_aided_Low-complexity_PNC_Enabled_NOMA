%  author: ph014@ie.cuhk.edu.hk
%  date:   April 20, 2015 / Modified 2020/11/16
%  demodulatoin

function [ ABits, BBits, CBits, DBits, X1Bits, X2Bits, X3Bits, X4Bits, X5Bits, X6Bits, X11Bits, X12Bits, X13Bits, X14Bits, ...
    X0Bits, p_diff_A,p_diff_B,p_diff_xor1 ] ...
    = NCMA_demod(rx, ha, hb, hc, hd, modulation, pkt_length, snr )

logp_A = zeros(pkt_length,1);
p_diff_A = zeros(pkt_length,1);
logp_B = zeros(pkt_length,1);
p_diff_B = zeros(pkt_length,1);
logp_C = zeros(pkt_length,1);
p_diff_C = zeros(pkt_length,1);
logp_D = zeros(pkt_length,1);
p_diff_D = zeros(pkt_length,1);

logp_xor1 = zeros(pkt_length,1);
p_diff_xor1 = zeros(pkt_length,1);
logp_xor2 = zeros(pkt_length,1);
p_diff_xor2 = zeros(pkt_length,1);
logp_xor3 = zeros(pkt_length,1);
p_diff_xor3 = zeros(pkt_length,1);
logp_xor4 = zeros(pkt_length,1);
p_diff_xor4 = zeros(pkt_length,1);
logp_xor5 = zeros(pkt_length,1);
p_diff_xor5 = zeros(pkt_length,1);
logp_xor6 = zeros(pkt_length,1);
p_diff_xor6 = zeros(pkt_length,1);
logp_xor11 = zeros(pkt_length,1);
p_diff_xor11 = zeros(pkt_length,1);
logp_xor12 = zeros(pkt_length,1);
p_diff_xor12 = zeros(pkt_length,1);
logp_xor13 = zeros(pkt_length,1);
p_diff_xor13 = zeros(pkt_length,1);
logp_xor14 = zeros(pkt_length,1);
p_diff_xor14 = zeros(pkt_length,1);
logp_xor0 = zeros(pkt_length,1);
p_diff_xor0 = zeros(pkt_length,1);

pA_0 = zeros(pkt_length,1);
pA_1 = zeros(pkt_length,1);
pB_0 = zeros(pkt_length,1);
pB_1 = zeros(pkt_length,1);
pC_0 = zeros(pkt_length,1);
pC_1 = zeros(pkt_length,1);
pD_0 = zeros(pkt_length,1);
pD_1 = zeros(pkt_length,1);
pX1_0 = zeros(pkt_length,1);
pX1_1 = zeros(pkt_length,1);
pX2_0 = zeros(pkt_length,1);
pX2_1 = zeros(pkt_length,1);
pX3_0 = zeros(pkt_length,1);
pX3_1 = zeros(pkt_length,1);
pX4_0 = zeros(pkt_length,1);
pX4_1 = zeros(pkt_length,1);
pX5_0 = zeros(pkt_length,1);
pX5_1 = zeros(pkt_length,1);
pX6_0 = zeros(pkt_length,1);
pX6_1 = zeros(pkt_length,1);
pX11_0 = zeros(pkt_length,1);
pX11_1 = zeros(pkt_length,1);
pX12_0 = zeros(pkt_length,1);
pX12_1 = zeros(pkt_length,1);
pX13_0 = zeros(pkt_length,1);
pX13_1 = zeros(pkt_length,1);
pX14_0 = zeros(pkt_length,1);
pX14_1 = zeros(pkt_length,1);
pX0_0 = zeros(pkt_length,1);
pX0_1 = zeros(pkt_length,1);

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
        d(pair_index) = norm( ha(ii) * data_pair(pair_index, 1)  ... 
            + hb(ii) * data_pair(pair_index,2)  ...
            + hc(ii) * data_pair(pair_index,3)  ...
            + hd(ii) * data_pair(pair_index,4)  ...
            - rx(ii) )^2;
    end
    
    if strcmp(modulation,'BPSK') == 1
        %packet A
        p0_A = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_A = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_A(ii) = p0_A-p1_A;
        logp_A(ii) = log(p0_A-p1_A);
        pA_0(ii) = p0_A;
        pA_1(ii) = p1_A;
        
        %packet B
        p0_B = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin)  ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_B = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_B(ii) = p0_B-p1_B;
        logp_B(ii) = log(p0_B-p1_B);
        pB_0(ii) = p0_B;
        pB_1(ii) = p1_B;
        
        %packet C
        p0_C = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_C = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin);
        p_diff_C(ii) = p0_C-p1_C;
        logp_C(ii) = log(p0_C-p1_C);
        pC_0(ii) = p0_C;
        pC_1(ii) = p1_C;
        
        %packet D
        p0_D = exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_D = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin);
        p_diff_D(ii) = p0_D-p1_D;
        logp_D(ii) = log(p0_D-p1_D);
        pD_0(ii) = p0_D;
        pD_1(ii) = p1_D;
        
        %packet xor1: AxorB
        p0_xor1 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor1 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor1(ii) = p0_xor1-p1_xor1;
        logp_xor1(ii) = log(p0_xor1-p1_xor1);
        pX1_0(ii) = p0_xor1;
        pX1_1(ii) = p1_xor1;
        
        %packet AxorC
        p0_xor2 = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor2 = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xor2(ii) = p0_xor2-p1_xor2;
        logp_xor2(ii) = log(p0_xor2-p1_xor2);
        pX2_0(ii) = p0_xor2;
        pX2_1(ii) = p1_xor2;
        
        %packet AxorD
        p0_xor3 = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor3 = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xor3(ii) = p0_xor3-p1_xor3;
        logp_xor3(ii) = log(p0_xor3-p1_xor3);
        pX3_0(ii) = p0_xor3;
        pX3_1(ii) = p1_xor3;
        
        %packet BxorC
        p0_xor4 = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor4 = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor4(ii) = p0_xor4-p1_xor4;
        logp_xor4(ii) = log(p0_xor4-p1_xor4);
        pX4_0(ii) = p0_xor4;
        pX4_1(ii) = p1_xor4;
        
        %packet BxorD
        p0_xor5 = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor5 = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor5(ii) = p0_xor5-p1_xor5;
        logp_xor5(ii) = log(p0_xor5-p1_xor5);
        pX5_0(ii) = p0_xor5;
        pX5_1(ii) = p1_xor5;
        
        %packet CxorD
        p0_xor6 = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin) ...
            +exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor6 = exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin);
        p_diff_xor6(ii) = p0_xor6-p1_xor6;
        logp_xor6(ii) = log(p0_xor6-p1_xor6);
        pX6_0(ii) = p0_xor6;
        pX6_1(ii) = p1_xor6;
        
        %packet AxorBxorC
        p0_xor11 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor11 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor11(ii) = p0_xor11-p1_xor11;
        logp_xor11(ii) = log(p0_xor11-p1_xor11);
        pX11_0(ii) = p0_xor11;
        pX11_1(ii) = p1_xor11;
        
        %packet AxorBxorD
        p0_xor12 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor12 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor12(ii) = p0_xor12-p1_xor12;
        logp_xor12(ii) = log(p0_xor12-p1_xor12);
        pX12_0(ii) = p0_xor12;
        pX12_1(ii) = p1_xor12;
        
        %packet AxorCxorD
        p0_xor13 = exp(-d(3)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(6)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(10)^2/snr_lin)+exp(-d(15)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor13 = exp(-d(1)^2/snr_lin)+exp(-d(2)^2/snr_lin)+exp(-d(7)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(11)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(14)^2/snr_lin);
        p_diff_xor13(ii) = p0_xor13-p1_xor13;
        logp_xor13(ii) = log(p0_xor13-p1_xor13);
        pX13_0(ii) = p0_xor13;
        pX13_1(ii) = p1_xor13;
        
        %packet BxorCxorD
        p0_xor14 = exp(-d(2)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor14 = exp(-d(1)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor14(ii) = p0_xor14-p1_xor14;
        logp_xor14(ii) = log(p0_xor14-p1_xor14);
        pX14_0(ii) = p0_xor14;
        pX14_1(ii) = p1_xor14;
        
        %packet AxorBxorCxorD
        p0_xor0 = exp(-d(1)^2/snr_lin)+exp(-d(4)^2/snr_lin)+exp(-d(6)^2/snr_lin)+exp(-d(7)^2/snr_lin) ...
            +exp(-d(10)^2/snr_lin)+exp(-d(11)^2/snr_lin)+exp(-d(13)^2/snr_lin)+exp(-d(16)^2/snr_lin);
        p1_xor0 = exp(-d(2)^2/snr_lin)+exp(-d(3)^2/snr_lin)+exp(-d(5)^2/snr_lin)+exp(-d(8)^2/snr_lin) ...
            +exp(-d(9)^2/snr_lin)+exp(-d(12)^2/snr_lin)+exp(-d(14)^2/snr_lin)+exp(-d(15)^2/snr_lin);
        p_diff_xor0(ii) = p0_xor0-p1_xor0;
        logp_xor0(ii) = log(p0_xor0-p1_xor0);
        pX0_0(ii) = p0_xor0;
        pX0_1(ii) = p1_xor0;
    end
    
    ABits = -1*reshape([pA_0,pA_1]',1,length(pA_0)*2); 
    BBits = -1*reshape([pB_0,pB_1]',1,length(pB_0)*2);
    CBits = -1*reshape([pC_0,pC_1]',1,length(pC_0)*2);
    DBits = -1*reshape([pD_0,pD_1]',1,length(pD_0)*2);
    X1Bits = -1*reshape([pX1_0,pX1_1]',1,length(pX1_0)*2);
    X2Bits = -1*reshape([pX2_0,pX2_1]',1,length(pX2_0)*2);
    X3Bits = -1*reshape([pX3_0,pX3_1]',1,length(pX3_0)*2);
    X4Bits = -1*reshape([pX4_0,pX4_1]',1,length(pX4_0)*2);
    X5Bits = -1*reshape([pX5_0,pX5_1]',1,length(pX5_0)*2);
    X6Bits = -1*reshape([pX6_0,pX6_1]',1,length(pX6_0)*2);
    X11Bits = -1*reshape([pX11_0,pX11_1]',1,length(pX11_0)*2);
    X12Bits = -1*reshape([pX12_0,pX12_1]',1,length(pX12_0)*2);
    X13Bits = -1*reshape([pX13_0,pX13_1]',1,length(pX13_0)*2);
    X14Bits = -1*reshape([pX14_0,pX14_1]',1,length(pX14_0)*2);
    X0Bits = -1*reshape([pX0_0,pX0_1]',1,length(pX0_0)*2);
end
end

