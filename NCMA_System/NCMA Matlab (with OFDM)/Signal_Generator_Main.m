
function [rxdata] = Signal_Generator_Main()

rxdata = [];

a = read_complex_binary('ncma_txdata1/txA_64S.dat');
b = read_complex_binary('ncma_txdata1/txB_64S.dat');

data = a(321:end);
Eb = norm(data)^2 /length(data);
L = length(a);% the pkt length
N = L * 3;

snr_a = 9.5;
sigma = sqrt(Eb/(10^(snr_a/10))); % noise variance
snr_d = 4;

b = b .* sqrt(power(10,snr_d/10));

rsnr = [];
s=[];
d=[];

cfo_a = 0; %100/1e6;
cfo_b = 0; %-100/1e6;

for i=1:1000
    r = (sigma/sqrt(2))*(randn(N,1) + 1j*randn(N,1));
    data_a = a(321:end);
    data_b = b(321:end);
    psnr_a=10*log10(norm(data_a)^2 / norm(r(1:length(data_a)))^2);
    psnr_b=10*log10(norm(data_b)^2 / norm(r(1:length(data_b)))^2);
    rsnr = [rsnr;psnr_a psnr_b];

    loc_s = 1*L + 1;
    
    num = randi([0, 360]);
    degree = num/360*2*pi;
    d=[d;degree];
    %r(loc_s : loc_s + L - 1) = r(loc_s:loc_s + L - 1) + a + b;
    %r(loc_s : loc_s + L - 1) = r(loc_s:loc_s + L - 1) + a + b * exp(1j*pi/4);
    r(loc_s : loc_s + L - 1) = r(loc_s:loc_s + L - 1) + a + b * exp(1j*degree);
    %r(loc_s : loc_s + L - 1) = a + b;

    %{
    % Rayleigh fading
    ha = (randn + 1j*randn)/sqrt(2);
    hb = (randn + 1j*randn)/sqrt(2);
    cfo_a_rotation = exp(1j*2*pi*cfo_a*[1:length(a)]');
    cfo_b_rotation = exp(1j*2*pi*cfo_b*[1:length(b)]');
    a_after_channel = ha * a .* cfo_a_rotation;
    b_after_channel = hb * b .* cfo_b_rotation;
    r(loc_s : loc_s + L - 1) = r(loc_s:loc_s + L - 1) + a_after_channel + b_after_channel;
    %}
    
    s=[s;r];
end
    
rsnr;
[mean(rsnr(:,1)) mean(rsnr(:,2))]
d(1:10:end)

%write_complex_binary(s,'awgnAB_0_64S_5dB.dat');
%write_complex_binary(s,'testAB_64S.dat');
%write_complex_binary(s,'awgnAB_random_64S_9dB.dat');
%write_complex_binary(s,'awgnAB_random_64S_75dB_135dB.dat');
%write_complex_binary(s,'awgnAB_random_64S_95dB_135dB.dat');
%write_complex_binary(s,'awgnAB_flat_64S_75dB_135dB.dat');
%write_complex_binary(s,'awgnAB_flat_64S_95dB_135dB.dat');
%x=read_complex_binary('testAB_64S.dat');
%Correlate(x,4);
%rxdata = s;


end
