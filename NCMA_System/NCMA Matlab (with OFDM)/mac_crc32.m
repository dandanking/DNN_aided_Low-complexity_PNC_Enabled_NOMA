% copyright: Daniel Halperin
% http://www.cs.washington.edu/homes/dhalperi/useful.html
%
% lzyou: it seems that crc32_ftw is different from crc32_halperin
% Take benchmark_data/pnc_pkt.16S.32f.1.dat for example
%   * cal_crc32_halperin=crc32(bits);
%   * cal_crc32_ftw=crc32_ftw(bits);
%
% Maybe crc_halperin is compatible with GNURadio, while crc_ftw not
% 

function ret = mac_crc32(bits)
poly = [1 de2bi(hex2dec('EDB88320'), 32)]';
bits = bits(:);

% Flip first 32 bits
bits(1:32) = 1 - bits(1:32);
% Add 32 zeros at the back
bits = [bits; zeros(32,1)];

% Initialize remainder to 0
rem = zeros(32,1);
% Main compution loop for the CRC32
for i = 1:length(bits)
    rem = [rem; bits(i)]; %#ok<AGROW>
    if rem(1) == 1
        rem = mod(rem + poly, 2);
    end
    rem = rem(2:33);
end

% Flip the remainder before returning it
ret = 1 - rem;
end
