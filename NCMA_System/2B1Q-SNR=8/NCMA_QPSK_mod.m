function [tx_ci,tx_cq] = NCMA_QPSK_mod(codeword_ci,codeword_cq)

%     I_c = codeword_c(1:2:length(codeword_c)-1);
%     Q_c = codeword_c(2:2:length(codeword_c));
%     tx_c = (2*I_c-1)+1j*(2*Q_c-1);
%     tx_c = tx_c/sqrt(2);
    tx_ci = 2*codeword_ci-1;
    tx_cq = 2*codeword_cq-1;
    
end

