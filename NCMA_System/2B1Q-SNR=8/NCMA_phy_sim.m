%  author: ph014@ie.cuhk.edu.hk
%          guanyu97@foxmail.com
%  date:   April 20, 2015 / Modified 2021/11/4
%  Main Program
%  Note: difference phase offsets between users affect decoding performance
   
clear

%matlabpool local 12
%trellis information
constraint_length = 7;%约束长度

trellis_encoded = poly2trellis(constraint_length,[133,171]);%卷积码生成多项式，编码规则/网格（2，1，7）
trellis_decoded = poly2trellis(constraint_length,[133,171]);

modulation = 'BPSK';

snr = 8:8;
degree_all = 45:45;

pkt_length = 64*48/2;
pkt_length_c = 64*48;
num_packet = 100000;

for de = 1:length(degree_all)
    for  pp = 1:length(snr)
        %degreeB = degree_all(de);
        fprintf('snr = %d\n', snr(pp));
        
        berA_all = zeros(num_packet,1);
        berB_all = zeros(num_packet,1);
        berC_all = zeros(num_packet,1);
        berX1_all = zeros(num_packet,1);
        berX2_all = zeros(num_packet,1);
        berX3_all = zeros(num_packet,1);
        berX0_all = zeros(num_packet,1);
        
        phy_raw_map = zeros(num_packet,13); % collect PHY-layer statistics
        
        for num = 1:num_packet
            %fprintf('packet = %d\n',num);
            degreeB = randi([0,359]);
            degreeC = randi([0,359]);
            msg_a = randi([0 1],pkt_length-8,1); % message sequence for A
            msg_a = [msg_a;zeros(8,1)];
            msg_b = randi([0 1],pkt_length-8,1); % message sequence for B
            msg_b = [msg_b;zeros(8,1)];
            msg_c = randi([0 1],pkt_length_c-16,1); % message sequence for C
            msg_c = [msg_c;zeros(16,1)];
            %拆分packet c
            msg_ci = msg_c(1:2:length(msg_c));
            msg_cq = msg_c(2:2:length(msg_c));
            
            msg_xorab = bitxor(msg_a,msg_b);
            msg_xoraci = bitxor(msg_a,msg_ci);
            msg_xoracq = bitxor(msg_a,msg_cq);
            msg_xorbci = bitxor(msg_b,msg_ci);
            msg_xorbcq = bitxor(msg_b,msg_cq);
            msg_xorabci = bitxor(msg_xorab,msg_ci);
            msg_xorabcq = bitxor(msg_xorab,msg_cq);
            
            codeword_a = convenc(msg_a,trellis_encoded);%卷积编码
            codeword_b = convenc(msg_b,trellis_encoded);
            codeword_ci = convenc(msg_ci,trellis_encoded);
            codeword_cq = convenc(msg_cq,trellis_encoded);
            
            %modulation
            [ tx_a, tx_b,~] = NCMA_mod(codeword_a, codeword_b, 0, modulation); %modulation=‘BPSK’，返回经调制后的码字
            [tx_ci,tx_cq] = NCMA_QPSK_mod(codeword_ci,codeword_cq);
            %channel
            n = 1/sqrt(2)*(randn(length(tx_a),1) + j*randn(length(tx_a),1));
            
            hci = ones(size(tx_ci));
            tx_hci = hci.*tx_ci;
            hcq = ones(size(tx_cq))*exp(j*pi*90/180);
            tx_hcq = hcq.*tx_cq;
            
            ha = ones(size(tx_a))*exp(j*pi*degreeB/180);
            tx_ha = ha.*tx_a;
            hb = ones(size(tx_b))*exp(j*pi*degreeC/180);             
            tx_hb = hb.*tx_b; 
                    
            % receiver starts
            rx = tx_hci+ tx_hcq + tx_ha + tx_hb + 10^(-snr(pp)/20)*n; 

            % demodulation 
            [logp_Ci, logp_Cq, logp_A, logp_B, logp_xorAB, logp_xorACi, logp_xorACq, logp_xorBCi, logp_xorBCq, logp_xorABCi, logp_xorABCq] ...
                = NCMA_demod(rx, hci, hcq, ha, hb, modulation, pkt_length, snr(pp));
               
            % decoder ber: convert to 0/1
            
             % encoded_A(find(encoded_A>=0))=0;
             % encoded_A(find(encoded_A<0))=1;
            
             % encoded_B(find(encoded_B>=0))=0;
             % encoded_B(find(encoded_B<0))=1;
            
             % encoded_X(find(encoded_X>=0))=0;
             % encoded_X(find(encoded_X<0))=1;
            
            decoded_source_A = viterbi_decoder_bit_level(logp_A,trellis_decoded,pkt_length);
            decoded_source_B = viterbi_decoder_bit_level(logp_B,trellis_decoded,pkt_length);
            decoded_source_Ci = viterbi_decoder_bit_level(logp_Ci,trellis_decoded,pkt_length);
            decoded_source_Cq = viterbi_decoder_bit_level(logp_Cq,trellis_decoded,pkt_length);
            decoded_source_xorAB = viterbi_decoder_bit_level(logp_xorAB,trellis_decoded,pkt_length);
            decoded_source_xorACi = viterbi_decoder_bit_level(logp_xorACi,trellis_decoded,pkt_length);
            decoded_source_xorACq = viterbi_decoder_bit_level(logp_xorACq,trellis_decoded,pkt_length);
            decoded_source_xorBCi = viterbi_decoder_bit_level(logp_xorBCi,trellis_decoded,pkt_length);
            decoded_source_xorBCq = viterbi_decoder_bit_level(logp_xorBCq,trellis_decoded,pkt_length);
            decoded_source_xorABCi = viterbi_decoder_bit_level(logp_xorABCi,trellis_decoded,pkt_length);
            decoded_source_xorABCq = viterbi_decoder_bit_level(logp_xorABCq,trellis_decoded,pkt_length);
            
            berA = sum(bitxor(decoded_source_A,msg_a))/length(msg_a);
            berB = sum(bitxor(decoded_source_B,msg_b))/length(msg_b);
            berCi = sum(bitxor(decoded_source_Ci,msg_ci))/length(msg_ci);
            berCq = sum(bitxor(decoded_source_Cq,msg_cq))/length(msg_cq);
            berXAB = sum(bitxor(decoded_source_xorAB,msg_xorab))/length(msg_xorab);
            berXACi = sum(bitxor(decoded_source_xorACi,msg_xoraci))/length(msg_xoraci);
            berXACq = sum(bitxor(decoded_source_xorACq,msg_xoracq))/length(msg_xoracq);
            berXBCi = sum(bitxor(decoded_source_xorBCi,msg_xorbci))/length(msg_xorbci);
            berXBCq = sum(bitxor(decoded_source_xorBCq,msg_xorbcq))/length(msg_xorbcq);
            berXABCi = sum(bitxor(decoded_source_xorABCi,msg_xorabci))/length(msg_xorabci);
            berXABCq = sum(bitxor(decoded_source_xorABCq,msg_xorabcq))/length(msg_xorabcq);
            
%             berA_all(num) = berA;
%             berB_all(num) = berB;
%             berCi_all(num) = berCi;
%             berCq_all(num) = berCq;
%             berX1_all(num) = berX1;
%             berX2_all(num) = berX2;
%             berX3_all(num) = berX3;
%             berX0_all(num) = berX0;
            
            okA=0; okB=0;okCi=0;okCq=0;okXAB=0;okXACi=0;okXACq=0;okXBCi=0;okXBCq=0;okXABCi=0;okXABCq=0;
            if berA == 0
                okA=1;
            end
            if berB == 0
                okB=1;
            end
            if berCi == 0
                okCi=1;
            end
            if berCq == 0
                okCq=1;
            end
            if berXAB == 0
                okXAB = 1;
            end
            if berXACi == 0
                okXACi = 1;
            end
            if berXACq == 0
                okXACq = 1;
            end
            if berXBCi == 0
                okXBCi = 1;
            end
            if berXBCq == 0
                okXBCq = 1;
            end
            if berXABCi == 0
                okXABCi = 1;
            end
            if berXABCq == 0
                okXABCq = 1;
            end
            phy_raw_map(num, :) = [degreeB, degreeC, okA, okB, okCi, okCq, okXAB, okXACi, okXACq, okXBCi, okXBCq, okXABCi, okXABCq ];
        end
        
    end
end

save('2B1Q-SNR=8.mat');
%matlabpool close
