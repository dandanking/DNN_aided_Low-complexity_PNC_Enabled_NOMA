
function [ok,ber,rxCRC,calCRC] = mac_crc32_wrapper(rx_bits,tx_bits,relay_flag)
%FIXME: should be aware of modulation and coding rate
%FIXME: add poly map automatically
poly_map = {'8CD5911F', '590CCDb6', 'A62A4EBA', '6EAF133A', '52FA50E2', 'EAB0E244'}; 


benchmark_nbits = 16*24;
number = ceil(length(rx_bits)/benchmark_nbits);
index = round(log(number)/log(2))+1;  % to make base log2
str=poly_map(index);

assert(length(str)>=0);
poly=de2bi(hex2dec(str), 32, 'left-msb')';

rx_raw_bits = rx_bits(1:length(rx_bits)-32);
cal_crc32=mac_crc32(rx_raw_bits);
rx_crc32 = rx_bits(length(rx_bits)-31:length(rx_bits));
if relay_flag
    rx_crc32 = xor(rx_crc32, poly);
end
if sum(bitxor(rx_crc32,cal_crc32))==0
    ok=1;
else
    ok=0;
end

rxCRC=rx_crc32;
calCRC=cal_crc32;

if relay_flag && sum(rx_raw_bits) == 0
    fprintf(' CRC WARNING: ALL ZEROS \n');
    ok=0;
end

ber=[];
if ~isempty(tx_bits)
    %tx_raw_bits=tx_bits(1:length(tx_bits)-32);
    %figure; stem(bitxor(rx_raw_bits,tx_raw_bits));
    ber=sum(bitxor(rx_bits,tx_bits))/length(tx_bits);
end