
function [RawBits, debug, debug2] = SIC_Data_Decoder(rxRawData, data_tones, nsym, H_a, H_b, knownTxRawData)
 % Return RawBits: [XOR, A, B]
 % Return debug: debug info per subcarrier

 counter = 1;
 debug=zeros(10000,10);
 debug2=zeros(10000,2);
 RawBits=zeros(data_tones*nsym,3);    
 
  
    for ii=1:nsym
        Data_Freq64 = rxRawData((ii-1)*64+1:ii*64);
        
        
        Pilot_phase1 = atan2( imag(Data_Freq64(33-21)./H_a(33-21)), real(Data_Freq64(33-21)./H_a(33-21)));
        Pilot_phase2 = atan2( imag(Data_Freq64(33-7)./H_b(33-7)),  real(Data_Freq64(33-7)./H_b(33-7)));
        Pilot_phase3 = atan2( imag(Data_Freq64(33+7)./H_a(33+7)),  real(Data_Freq64(33+7)./H_a(33+7)));
        Pilot_phase4 = atan2( imag(Data_Freq64(33+21)./H_b(33+21)), real(Data_Freq64(33+21)./H_b(33+21)));
        
        rotation_a = exp(1j*((Pilot_phase1+Pilot_phase3)/2));
        rotation_b = exp(1j*((Pilot_phase2+Pilot_phase4)/2));
        
        Data_buffer=zeros(64,3);
        
        range = [7:11 13:25 27:32 34:39 41:53 55:59];
        for d = range;
            HA = H_a(d); HB = H_b(d); in = Data_Freq64(d);
            %debug(counter,:) = [real(in) imag(in) real(HA) imag(HA) real(rotation_a) imag(rotation_a) real(HB) imag(HB) real(rotation_b) imag(rotation_b)];            
            HA = HA  * rotation_a;
            HB = HB  * rotation_b;
            % debug2(counter,:)
            [Data_buffer(d,:),~] = SIC_Subcarrier_Demod(in,HA,HB,H_a,H_b,knownTxRawData(counter,:));
            counter = counter+1;
        end
        
        RawBits((ii-1)*48+1 : ii*48, : ) = Data_buffer(range,:);        
        
        for d=1:64
            H_a(d) = H_a(d) * rotation_a;
            H_b(d) = H_b(d) * rotation_b;
        end
    end
    debug=debug(1:counter-1,:);
