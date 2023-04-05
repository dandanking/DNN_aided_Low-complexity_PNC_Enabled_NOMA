


close all;
clear all;

% ------------------------------------------------------------------------%
%                         Common Parameters                               %
% ------------------------------------------------------------------------%
global DECODER_TYPE
fft_size=64;
data_tones=48;
padding_bytes=1;
nsym=64;
nbits=nsym*data_tones;

Num_data_symbols=nsym;

DECODER_TYPE = 'soft';
mode = 3; %FIXME: change the code to be workable when mode=3

% ------------------------------------------------------------------------%
%                          Tx Data Part                                   %
% ------------------------------------------------------------------------%
strA = sprintf('ncma_txdata2/txA_%dS_intrlv.datb',Num_data_symbols);
strB = sprintf('ncma_txdata2/txB_%dS_intrlv.datb',Num_data_symbols);
txA_encoded_bits=read_char_binary(strA);
txB_encoded_bits=read_char_binary(strB);
txX_encoded_bits=bitxor(txA_encoded_bits,txB_encoded_bits);
txA_source_bits=phy_viterbi_decoder(phy_de_interleaver(txA_encoded_bits),'hard');
txB_source_bits=phy_viterbi_decoder(phy_de_interleaver(txB_encoded_bits),'hard');
txX_source_bits=phy_viterbi_decoder(phy_de_interleaver(txX_encoded_bits),'hard');
% remove padding bits
txA_source_bits=txA_source_bits(1:end-8);
txB_source_bits=txB_source_bits(1:end-8);
txX_source_bits=txX_source_bits(1:end-8);

% Mode 3: begin from soft demod output
% load GNURadio's soft output
%%{
rxX=read_uint8_binary('debug_data/ncma-x-bits.dat');
rxA=read_uint8_binary('debug_data/ncma-a-bits.dat');
rxB=read_uint8_binary('debug_data/ncma-b-bits.dat');
%}

for jj = 1:1
    rawXBits = rxX((jj-1)*nbits+1:jj*nbits);
    rawABits = rxA((jj-1)*nbits+1:jj*nbits);
    rawBBits = rxB(((jj-1)*nbits+1:jj*nbits));

    XBits = phy_viterbi_decoder(phy_de_interleaver(rawXBits),DECODER_TYPE);
    ABits = phy_viterbi_decoder(phy_de_interleaver(rawABits),DECODER_TYPE);
    BBits = phy_viterbi_decoder(phy_de_interleaver(rawBBits),DECODER_TYPE);
    % remove padding bits
    XBits=XBits(1:end-8);
    ABits=ABits(1:end-8);
    BBits=BBits(1:end-8);

end