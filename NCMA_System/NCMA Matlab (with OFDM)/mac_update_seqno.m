%% Given count seq information, do update
%  we copy codes from previous version
%  hence we preserve all the names
%  author: lzyou@inc.cuhk.edu.hk
%  date:   April 12, 2013

function [outputSeqInfo,outputMap,d_X_decodable] = mac_update_seqno(cur_index,inputSeqInfo,d_X_decodable,mapPktIndex_ncma,L_A,L_B,ncma_mode, ...
                                                                    rxRawData,data_tones,nsym,H_a,H_b, ...
                                                                    txA_encoded_bits,txB_encoded_bits,txX_encoded_bits, ...
                                                                    txA_source_bits,txB_source_bits,txX_source_bits)
    % extract info from input parameters
    A_INDEX = 1; B_INDEX = 2; X_INDEX = 3;
    countA_ncma = inputSeqInfo(A_INDEX,1); pktAnum_ncma = inputSeqInfo(A_INDEX,2);
    countB_ncma = inputSeqInfo(B_INDEX,1); pktBnum_ncma = inputSeqInfo(B_INDEX,2);
    countX_ncma = inputSeqInfo(X_INDEX,1); pktXnum_ncma = inputSeqInfo(X_INDEX,2);
    
    if strcmp(ncma_mode,'ncma5') || strcmp(ncma_mode,'ncma6')
        mode = 'back';
    else
        mode = [];
    end
    
    global fphy
    global fmac
          
    % common definition
    TRY_XOR_EQ = 1; % in this version, we need *TRY XOR*
    NULL_SIGNAL = 0;
    SELF_RECOVER_SIGNAL = -1;
    XOR_RECOVER_SIGNAL  = -10;  % must be consist with main function
    L_X = max(L_A,L_B);
    
    null_encoded_bits=-10*ones(length(txA_encoded_bits),1);
    
    % fill in XOR recovery information
    if ~mapPktIndex_ncma(cur_index,X_INDEX) && d_X_decodable
        mapPktIndex_ncma(cur_index,X_INDEX) = SELF_RECOVER_SIGNAL * pktXnum_ncma;
    end
    
    % assumption: either A -> B, or B -> A. Iterative is impossible.
    % while is used for update
    d_updateA_indicator = 0;
    d_updateB_indicator = 0;
    d_count = 0;        
    while ((countA_ncma >= L_A) || (countB_ncma >= L_B) || (TRY_XOR_EQ && (d_X_decodable == 0) && (countX_ncma >= L_X)))
        %fprintf(fmac,'%s:  d_count = %d curIndex=%d countA=[%d %d] countB=[%d %d] countX=[%d %d] \n',ncma_mode,d_count,cur_index,pktAnum_ncma,countA_ncma,pktBnum_ncma,countB_ncma,pktXnum_ncma,countX_ncma);
        last_index = max(find( mapPktIndex_ncma(:,A_INDEX) == (pktAnum_ncma-1) | ...
                               mapPktIndex_ncma(:,B_INDEX) == (pktBnum_ncma-1) | ...
                               mapPktIndex_ncma(:,X_INDEX) == (pktXnum_ncma-1) | ...
                               mapPktIndex_ncma(:,A_INDEX) == SELF_RECOVER_SIGNAL*(pktAnum_ncma-1) | ...
                               mapPktIndex_ncma(:,B_INDEX) == SELF_RECOVER_SIGNAL*(pktBnum_ncma-1) | ...
                               mapPktIndex_ncma(:,X_INDEX) == SELF_RECOVER_SIGNAL*(pktXnum_ncma-1) | ...
                               mapPktIndex_ncma(:,A_INDEX) == XOR_RECOVER_SIGNAL *(pktAnum_ncma-1) | ...
                               mapPktIndex_ncma(:,B_INDEX) == XOR_RECOVER_SIGNAL *(pktBnum_ncma-1) | ...
                               mapPktIndex_ncma(:,X_INDEX) == XOR_RECOVER_SIGNAL *(pktXnum_ncma-1)));
        
        d_count = d_count + 1;
        assert(d_count <= 3);                               % it is possible that X -> B -> A
        
        if countA_ncma >= L_A
            d_updateA_indicator = 1;                        % we have updated A
            
            if d_updateB_indicator == 0
                if pktAnum_ncma == 1
                    last_index = 0;
                end
                for pkt_index = last_index+1:cur_index
                    if mapPktIndex_ncma(pkt_index, A_INDEX) == NULL_SIGNAL
                        mapPktIndex_ncma(pkt_index, A_INDEX) = SELF_RECOVER_SIGNAL * pktAnum_ncma;
                    end
                    
                    % only have XOR packet, but now we recover A, update B
                    % two cases: got X, or X is recoverable
                    if (mapPktIndex_ncma(pkt_index, X_INDEX) || (TRY_XOR_EQ && d_X_decodable)) && (~mapPktIndex_ncma(pkt_index,B_INDEX))
                        countB_ncma = countB_ncma + 1;
                        mapPktIndex_ncma(pkt_index,B_INDEX) = XOR_RECOVER_SIGNAL * pktBnum_ncma;
                    end
                    
                    if strcmp(mode,'back')==1 && ~mapPktIndex_ncma(pkt_index, B_INDEX) && ~mapPktIndex_ncma(pkt_index, X_INDEX) && countB_ncma < L_B
                        % Nothing is received. Given A, we try to recover
                        % B using two-phase decoding.
                        knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
                        knownTxRawData(:,A_INDEX) = txA_encoded_bits;
                        RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                        rawXBits = RawBits(:,X_INDEX); rawBBits = RawBits(:,B_INDEX);
                        XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),'soft');
                        BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),'soft');
                        [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits(1:end-8),txX_source_bits,1);
                        [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
                        if okX || okB
                            countB_ncma = countB_ncma + 1;
                            mapPktIndex_ncma(pkt_index,B_INDEX) = pktBnum_ncma; %XOR_RECOVER_SIGNAL * pktBnum_ncma;
                            countX_ncma = countX_ncma + 1;
                            mapPktIndex_ncma(pkt_index,X_INDEX) = pktXnum_ncma; %XOR_RECOVER_SIGNAL * pktXnum_ncma;
                        end
                        fprintf(fphy,'%s:  backward A+%d-%d: okX=%d berX=%f | okB=%d berB=%f | countA=[%d %d] countB=[%d %d] \n', ncma_mode, pkt_index, cur_index, okX, berX, okB, berB, pktAnum_ncma, countA_ncma, pktBnum_ncma, countB_ncma);
                    end
                end
            end % end of d_updateB_indicator
            
            pktAnum_ncma = pktAnum_ncma + 1;
            countA_ncma = 0;
            
            % if B has not been updated (ie: A -> B), also update XOR
            if d_updateB_indicator == 0
                d_X_decodable = 0;
                pktXnum_ncma = pktXnum_ncma + 1;
                countX_ncma = 0;
            end
        end
        
        if countB_ncma >= L_B
            d_updateB_indicator = 1;                        % we have updated B
            
            if d_updateA_indicator == 0     
                if pktBnum_ncma == 1
                    last_index = 0;
                end
                for pkt_index = last_index+1:cur_index
                    if mapPktIndex_ncma(pkt_index, B_INDEX) == NULL_SIGNAL
                        mapPktIndex_ncma(pkt_index, B_INDEX) = SELF_RECOVER_SIGNAL * pktBnum_ncma;
                    end
                    
                    % only have XOR packet, but now we recover A, update B
                    % two cases: got X, or X is recoverable
                    if pkt_index == 60
                        pkt_index = 60;
                    end
                    
                    % recover A when A is invalid and X is valid, no matter
                    % B is valid or not
                    if ((TRY_XOR_EQ && d_X_decodable) || mapPktIndex_ncma(pkt_index, X_INDEX)) && (~mapPktIndex_ncma(pkt_index,A_INDEX))
                        countA_ncma = countA_ncma + 1;
                        mapPktIndex_ncma(pkt_index,A_INDEX) = XOR_RECOVER_SIGNAL * pktAnum_ncma;
                    end
                    
                    % try two-phase decoding: 1) backward mode; 2) A is
                    % invalid; 3) B can be recovered right now.
                    if strcmp(mode,'back')==1 && ~mapPktIndex_ncma(pkt_index, A_INDEX) && ~mapPktIndex_ncma(pkt_index, X_INDEX) && countA_ncma < L_A
                        knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
                        knownTxRawData(:,B_INDEX) = txB_encoded_bits;
                        RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                        rawXBits = RawBits(:,X_INDEX); rawABits = RawBits(:,A_INDEX);
                        XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),'soft');
                        ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),'soft');
                        [okX,berX,rx_crc32,cal_crc32] = mac_crc32_wrapper(XBits(1:end-8),txX_source_bits,1);
                        [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
                        if okX || okA
                            countA_ncma = countA_ncma + 1;
                            mapPktIndex_ncma(pkt_index,A_INDEX) = pktAnum_ncma; %XOR_RECOVER_SIGNAL * pktAnum_ncma;
                            countX_ncma = countX_ncma + 1;
                            mapPktIndex_ncma(pkt_index,X_INDEX) = pktXnum_ncma; %XOR_RECOVER_SIGNAL * pktXnum_ncma;
                        end
                        fprintf(fphy,'%s:  backward B+%d-%d: okX=%d berX=%f | okA=%d berA=%f | countA=[%d %d] countB=[%d %d]  \n', ncma_mode, pkt_index, cur_index, okX, berX, okA, berA, pktAnum_ncma, countA_ncma, pktBnum_ncma, countB_ncma);
                    end
                end
            end % end of d_updateA_indicator
            
            pktBnum_ncma = pktBnum_ncma + 1;
            countB_ncma = 0;
            
            % if A has not been updated (ie: B -> A), also update XOR
            if d_updateA_indicator == 0                    
                d_X_decodable = 0;
                pktXnum_ncma = pktXnum_ncma + 1;
                countX_ncma = 0;
            end
        end % end of B
        
        if TRY_XOR_EQ && d_X_decodable == 0 && countX_ncma >= L_X
            d_X_decodable = 1;
            
            if d_updateA_indicator == 0 && d_updateB_indicator == 0
                if pktXnum_ncma == 1
                    last_index = 0;
                end
                
                % Recover X, and then A,B
                for pkt_index = last_index+1:cur_index
                    if mapPktIndex_ncma(pkt_index, X_INDEX) == NULL_SIGNAL
                        mapPktIndex_ncma(pkt_index, X_INDEX) = SELF_RECOVER_SIGNAL * pktXnum_ncma;
                        
                        if mapPktIndex_ncma(pkt_index,B_INDEX) && ~mapPktIndex_ncma(pkt_index,A_INDEX)
                            countA_ncma = countA_ncma + 1;
                            mapPktIndex_ncma(pkt_index,A_INDEX) = XOR_RECOVER_SIGNAL * pktAnum_ncma; %
                            %equations_x_gain = equations_x_gain + 1;
                        elseif mapPktIndex_ncma(pkt_index,A_INDEX) && ~mapPktIndex_ncma(pkt_index,B_INDEX)
                            countB_ncma = countB_ncma + 1;
                            mapPktIndex_ncma(pkt_index,B_INDEX) = XOR_RECOVER_SIGNAL * pktBnum_ncma; %
                            %equations_x_gain = equations_x_gain + 1;
                        elseif strcmp(mode,'back')==1 && ~mapPktIndex_ncma(pkt_index,A_INDEX) && ~mapPktIndex_ncma(pkt_index,B_INDEX) && countA_ncma < L_A && countB_ncma < L_B
                            % two-phase deoding
                            knownTxRawData = [null_encoded_bits,null_encoded_bits,null_encoded_bits];
                            knownTxRawData(:,X_INDEX) = txX_encoded_bits;
                            RawBits = PNC_Data_Decoder(rxRawData,data_tones,nsym,H_a,H_b,knownTxRawData);
                            rawABits = RawBits(:,A_INDEX); rawBBits = RawBits(:,B_INDEX);
                            ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),'soft');
                            BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),'soft');
                            [okA,berA,rx_crc32,cal_crc32] = mac_crc32_wrapper(ABits(1:end-8),txA_source_bits,0);
                            [okB,berB,rx_crc32,cal_crc32] = mac_crc32_wrapper(BBits(1:end-8),txB_source_bits,0);
                            if okA || okB
                                countA_ncma = countA_ncma + 1;
                                mapPktIndex_ncma(pkt_index,A_INDEX) = pktAnum_ncma; %XOR_RECOVER_SIGNAL * pktAnum_ncma;
                                countB_ncma = countB_ncma + 1;
                                mapPktIndex_ncma(pkt_index,B_INDEX) = pktBnum_ncma; %XOR_RECOVER_SIGNAL * pktBnum_ncma;
                            end
                            fprintf(fphy,'%s:  backward X+%d-%d: okA=%d berA=%f | okB=%d berB=%f | countA=[%d %d] countB=[%d %d] \n', ncma_mode, pkt_index, cur_index, okA, berA, okB, berB, pktAnum_ncma, countA_ncma, pktBnum_ncma, countB_ncma);
                        end
                    end
                end
            end
        end % end of X
    end % end of while

    % wrap info for output parameters
    outputMap = mapPktIndex_ncma;
    outputSeqInfo = zeros(3,2);
    outputSeqInfo(A_INDEX,:) = [countA_ncma, pktAnum_ncma];
    outputSeqInfo(B_INDEX,:) = [countB_ncma, pktBnum_ncma];
    outputSeqInfo(X_INDEX,:) = [countX_ncma, pktXnum_ncma];
    
end % end of function