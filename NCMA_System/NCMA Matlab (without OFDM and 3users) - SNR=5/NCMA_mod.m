%  author: ph014@ie.cuhk.edu.hk
%  date:   April 20, 2015 / Modified 2020/11/16
%  Modulation: gray mapping
%  Match GNURadio Implementation 

function [ tx_a, tx_b, tx_c ] = NCMA_mod(codeword_a, codeword_b, codeword_c, modulation) %function [输出变量] = 函数名称(输入变量）

if strcmp(modulation,'BPSK') == 1
    tx_a = 2*codeword_a-1;  %0 map to -1; 1 map to 1
    tx_b = 2*codeword_b-1;
    tx_c = 2*codeword_c-1;
    
elseif strcmp(modulation,'QPSK') == 1
    I_a = codeword_a(1:2:length(codeword_a)-1);
    Q_a = codeword_a(2:2:length(codeword_a));
    tx_a = (2*I_a-1)+1j*(2*Q_a-1);
    tx_a = tx_a/sqrt(2);
    
    I_b = codeword_b(1:2:length(codeword_b)-1);
    Q_b = codeword_b(2:2:length(codeword_b));
    tx_b = (2*I_b-1)+1j*(2*Q_b-1);
    tx_b = tx_b/sqrt(2);
end

end

