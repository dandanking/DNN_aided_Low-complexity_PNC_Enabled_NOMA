%  author: ph014@ie.cuhk.edu.hk
%          guanyu97@foxmail.com
%  date:   April 20, 2015 / Modified 2020/11/16
%  Main Program
%  Note: difference phase offsets between users affect decoding performance
   
clear

%matlabpool local 12
%trellis information
constraint_length = 7;%约束长度

trellis_encoded = poly2trellis(constraint_length,[133,171]);%卷积码生成多项式，编码规则/网格（2，1，7）
trellis_decoded = poly2trellis(constraint_length,[133,171]);

modulation = 'BPSK';
%modulation = 'QPSK';

snr = 8:8;
degree_all = 45:45;

pkt_length = 64*48/2;
num_packet = 100000;

for de = 1:length(degree_all)
    for  pp = 1:length(snr)
        %degreeB = degree_all(de);
        fprintf('snr = %d\n', snr(pp));
        
        berA_all = zeros(num_packet,1);
        berB_all = zeros(num_packet,1);
        berC_all = zeros(num_packet,1);
        berD_all = zeros(num_packet,1);
        berX1_all = zeros(num_packet,1);
        berX2_all = zeros(num_packet,1);
        berX3_all = zeros(num_packet,1);
        berX0_all = zeros(num_packet,1);

        
        phy_raw_map = zeros(num_packet,18); % collect PHY-layer statistics
        
        for num = 1:num_packet
            %fprintf('packet = %d\n',num);
            degreeB = randi([0,360]);
            degreeC = randi([0,360]);
            degreeD = randi([0,360]);
            msg_a = randi([0 1],pkt_length-8,1); % message sequence for A
            msg_a = [msg_a;zeros(8,1)];
            msg_b = randi([0 1],pkt_length-8,1); % message sequence for B
            msg_b = [msg_b;zeros(8,1)];
            msg_c = randi([0 1],pkt_length-8,1); % message sequence for C
            msg_c = [msg_c;zeros(8,1)];
            msg_d = randi([0 1],pkt_length-8,1); % message sequence for D
            msg_d = [msg_d;zeros(8,1)];
            msg_xor1 = bitxor(msg_a,msg_b);
            msg_xor2 = bitxor(msg_a,msg_c);
            msg_xor3 = bitxor(msg_a,msg_d);
            msg_xor4 = bitxor(msg_b,msg_c);
            msg_xor5 = bitxor(msg_b,msg_d);
            msg_xor6 = bitxor(msg_c,msg_d);
            msg_xor11 = bitxor(msg_xor1,msg_c);
            msg_xor12 = bitxor(msg_xor1,msg_d);
            msg_xor13 = bitxor(msg_xor2,msg_d);
            msg_xor14 = bitxor(msg_xor4,msg_d);
            msg_xor0 = bitxor(msg_xor11,msg_d);
            
            codeword_a = convenc(msg_a,trellis_encoded);%卷积编码
            codeword_b = convenc(msg_b,trellis_encoded);
            codeword_c = convenc(msg_c,trellis_encoded);
            codeword_d = convenc(msg_d,trellis_encoded);
            
            %modulation
            [ tx_a, tx_b, tx_c, tx_d] = NCMA_mod(codeword_a, codeword_b, codeword_c, codeword_d, modulation); %modulation=‘BPSK’，返回经调制后的码字
            
            %channel
            n = 1/sqrt(2)*(randn(length(tx_a),1) + j*randn(length(tx_a),1)); %返回一个均值为0，方差为1/2+j的length(tx_a)*1的随机数矩阵
            
            ha=ones(size(tx_a));
            tx_ha = ha.*tx_a;
                        
            hb=ones(size(tx_b))*exp(j*pi*degreeB/180);             
            tx_hb = hb.*tx_b; 
            
            hc=ones(size(tx_c))*exp(j*pi*degreeC/180);
            tx_hc = hc.*tx_c;
            
            hd=ones(size(tx_d))*exp(j*pi*degreeD/180);
            tx_hd = hd.*tx_d;
            
            % receiver starts
            rx = tx_ha + tx_hb+ tx_hc + tx_hd+ 10^(-snr(pp)/20)*n; 

            
            %demodulation 
            [ logp_A, logp_B, logp_C, logp_D, logp_xor1, logp_xor2, logp_xor3, logp_xor4, logp_xor5, logp_xor6, ... 
                logp_xor11, logp_xor12, logp_xor13, logp_xor14, logp_xor0, ... 
                encoded_A, encoded_B, encoded_X ] = NCMA_demod(rx, ha, hb, hc, hd, modulation, pkt_length, snr(pp));
               
            %decoder ber: convert to 0/1
            
             % encoded_A(find(encoded_A>=0))=0;
             % encoded_A(find(encoded_A<0))=1;
            
             % encoded_B(find(encoded_B>=0))=0;
             % encoded_B(find(encoded_B<0))=1;
            
             % encoded_X(find(encoded_X>=0))=0;
             % encoded_X(find(encoded_X<0))=1;
            
            decoded_source_A = viterbi_decoder_bit_level(logp_A,trellis_decoded,pkt_length);
            decoded_source_B = viterbi_decoder_bit_level(logp_B,trellis_decoded,pkt_length);
            decoded_source_C = viterbi_decoder_bit_level(logp_C,trellis_decoded,pkt_length);
            decoded_source_D = viterbi_decoder_bit_level(logp_D,trellis_decoded,pkt_length);
            decoded_source_xor1 = viterbi_decoder_bit_level(logp_xor1,trellis_decoded,pkt_length);
            decoded_source_xor2 = viterbi_decoder_bit_level(logp_xor2,trellis_decoded,pkt_length);
            decoded_source_xor3 = viterbi_decoder_bit_level(logp_xor3,trellis_decoded,pkt_length);
            decoded_source_xor4 = viterbi_decoder_bit_level(logp_xor4,trellis_decoded,pkt_length);
            decoded_source_xor5 = viterbi_decoder_bit_level(logp_xor5,trellis_decoded,pkt_length);
            decoded_source_xor6 = viterbi_decoder_bit_level(logp_xor6,trellis_decoded,pkt_length);
            decoded_source_xor11 = viterbi_decoder_bit_level(logp_xor11,trellis_decoded,pkt_length);
            decoded_source_xor12 = viterbi_decoder_bit_level(logp_xor12,trellis_decoded,pkt_length);
            decoded_source_xor13 = viterbi_decoder_bit_level(logp_xor13,trellis_decoded,pkt_length);
            decoded_source_xor14 = viterbi_decoder_bit_level(logp_xor14,trellis_decoded,pkt_length);
            decoded_source_xor0 = viterbi_decoder_bit_level(logp_xor0,trellis_decoded,pkt_length);
            
            berA = sum(bitxor(decoded_source_A,msg_a))/length(msg_a);
            berB = sum(bitxor(decoded_source_B,msg_b))/length(msg_b);
            berC = sum(bitxor(decoded_source_C,msg_c))/length(msg_c);
            berD = sum(bitxor(decoded_source_D,msg_d))/length(msg_d);
            berX1 = sum(bitxor(decoded_source_xor1,msg_xor1))/length(msg_xor1);
            berX2 = sum(bitxor(decoded_source_xor2,msg_xor2))/length(msg_xor2);
            berX3 = sum(bitxor(decoded_source_xor3,msg_xor3))/length(msg_xor3);
            berX4 = sum(bitxor(decoded_source_xor4,msg_xor4))/length(msg_xor4);
            berX5 = sum(bitxor(decoded_source_xor5,msg_xor5))/length(msg_xor5);
            berX6 = sum(bitxor(decoded_source_xor6,msg_xor6))/length(msg_xor6);
            berX11 = sum(bitxor(decoded_source_xor11,msg_xor11))/length(msg_xor11);
            berX12 = sum(bitxor(decoded_source_xor12,msg_xor12))/length(msg_xor12);
            berX13 = sum(bitxor(decoded_source_xor13,msg_xor13))/length(msg_xor13);
            berX14 = sum(bitxor(decoded_source_xor14,msg_xor14))/length(msg_xor14);
            berX0 = sum(bitxor(decoded_source_xor0,msg_xor0))/length(msg_xor0);
            
            berA_all(num) = berA;
            berB_all(num) = berB;
            berC_all(num) = berC;
            berX1_all(num) = berX1;
            berX2_all(num) = berX2;
            berX3_all(num) = berX3;
            berX0_all(num) = berX0;
            
            okA=0;okB=0;okC=0;okD=0;okX1=0;okX2=0;okX3=0;okX4=0;okX5=0;okX6=0;okX11=0;okX12=0;okX13=0;okX14=0;okX0=0;
            if berA == 0
                okA=1;
            end
            if berB == 0
                okB=1;
            end
            if berC == 0
                okC=1;
            end
            if berD == 0
                okD=1;
            end
            if berX1 == 0
                okX1 = 1;
            end
            if berX2 == 0
                okX2 = 1;
            end
            if berX3 == 0
                okX3 = 1;
            end
            if berX4 == 0
                okX4 = 1;
            end
            if berX5 == 0
                okX5 = 1;
            end
            if berX6 == 0
                okX6 = 1;
            end
            if berX11 == 0
                okX11 = 1;
            end
            if berX12 == 0
                okX12 = 1;
            end
            if berX13 == 0
                okX13 = 1;
            end
            if berX14 == 0
                okX14 = 1;
            end
            if berX0 == 0
                okX0 = 1;
            end
            
            phy_raw_map(num, :) = [degreeB, degreeC, degreeD, okA, okB, okC, okD, okX1, okX2, okX3, okX4, okX5, okX6, okX11, okX12,  okX13, okX14, okX0];
        end
        
    end
end


%matlabpool close
