%  author: ph014@ie.cuhk.edu.hk
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

snr = 5:5;
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
        berX1_all = zeros(num_packet,1);
        berX2_all = zeros(num_packet,1);
        berX3_all = zeros(num_packet,1);
        berX0_all = zeros(num_packet,1);

        
        phy_raw_map = zeros(num_packet,9); % collect PHY-layer statistics
        
        for num = 1:num_packet
            %fprintf('packet = %d\n',num);
            degreeB = randi([0,360]);
            degreeC = randi([0,360]);
            msg_a = randi([0 1],pkt_length-8,1); % message sequence for A
            msg_a = [msg_a;zeros(8,1)];
            msg_b = randi([0 1],pkt_length-8,1); % message sequence for B
            msg_b = [msg_b;zeros(8,1)];
            msg_c = randi([0 1],pkt_length-8,1); % message sequence for A
            msg_c = [msg_c;zeros(8,1)];
            msg_xor1 = bitxor(msg_a,msg_b);
            msg_xor2 = bitxor(msg_a,msg_c);
            msg_xor3 = bitxor(msg_b,msg_c);
            msg_xor0 = bitxor(msg_xor1,msg_c);
            
            codeword_a = convenc(msg_a,trellis_encoded);%卷积编码
            codeword_b = convenc(msg_b,trellis_encoded);
            codeword_c = convenc(msg_c,trellis_encoded);
            
            %modulation
            [ tx_a, tx_b, tx_c] = NCMA_mod(codeword_a, codeword_b, codeword_c, modulation); %modulation=‘BPSK’，返回经调制后的码字
            
            %channel
            n = 1/sqrt(2)*(randn(length(tx_a),1) + j*randn(length(tx_a),1)); %返回一个均值为0，方差为1/2+j的length(tx_a)*1的随机数矩阵
            
            ha=ones(size(tx_a));
            tx_ha = ha.*tx_a;
            
            hc=ones(size(tx_c))*exp(j*pi*degreeC/180);
            tx_hc = hc.*tx_c;
                             
            hb=ones(size(tx_b))*exp(j*pi*degreeB/180);             
            tx_hb = hb.*tx_b; 
                                
            % receiver starts
            rx = tx_ha + tx_hb+ tx_hc + 10^(-snr(pp)/20)*n; 

            
            %demodulation 
            [ logp_A, logp_B, logp_C, logp_xor1, logp_xor2, logp_xor3, logp_xor0, encoded_A, encoded_B, encoded_X ] = NCMA_demod(rx, ha, hb, hc, modulation, pkt_length, snr(pp));
               
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
            decoded_source_xor1 = viterbi_decoder_bit_level(logp_xor1,trellis_decoded,pkt_length);
            decoded_source_xor2 = viterbi_decoder_bit_level(logp_xor2,trellis_decoded,pkt_length);
            decoded_source_xor3 = viterbi_decoder_bit_level(logp_xor3,trellis_decoded,pkt_length);
            decoded_source_xor0 = viterbi_decoder_bit_level(logp_xor0,trellis_decoded,pkt_length);
            
            berA = sum(bitxor(decoded_source_A,msg_a))/length(msg_a);
            berB = sum(bitxor(decoded_source_B,msg_b))/length(msg_b);
            berC = sum(bitxor(decoded_source_C,msg_c))/length(msg_c);
            berX1 = sum(bitxor(decoded_source_xor1,msg_xor1))/length(msg_xor1);
            berX2 = sum(bitxor(decoded_source_xor2,msg_xor2))/length(msg_xor2);
            berX3 = sum(bitxor(decoded_source_xor3,msg_xor3))/length(msg_xor3);
            berX0 = sum(bitxor(decoded_source_xor0,msg_xor0))/length(msg_xor0);
            
            berA_all(num) = berA;
            berB_all(num) = berB;
            berC_all(num) = berC;
            berX1_all(num) = berX1;
            berX2_all(num) = berX2;
            berX3_all(num) = berX3;
            berX0_all(num) = berX0;
            
            okA=0; okB=0;okC=0;okX1=0;okX2=0;okX3=0;okX0=0;
            if berA == 0
                okA=1;
            end
            if berB == 0
                okB=1;
            end
            if berC == 0
                okC=1;
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
            if berX0 == 0
                okX0 = 1;
            end
            
            phy_raw_map(num, :) = [degreeB, degreeC, okA, okB, okC, okX1, okX2, okX3, okX0];
        end
        
    end
end


%matlabpool close
