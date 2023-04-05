
function output_pair = PNC_Subcarrier_Demod(in,H_A,H_B,H_A_list,H_B_list, knownTxRawData)
% Input:
%   in: recived point
%   H_A/H_B: channel informaion of user A/B
%   H_A_list/H_B_list: channel of all subcarriers
%   knownTxRawData: tx bits, for cross-layer decoding
% Output:
%   output: decoding bits (soft or hard), order=[A, B, X]
% TODO:
%   fix SIC algorithms for QPSK (July 29, 2014)

global DECODER_TYPE

VALID_INDEX = [7:32 34:59];
maxHA = max(abs(H_A_list(VALID_INDEX)));
maxHB = max(abs(H_B_list(VALID_INDEX)));
maxHAB = 0;
for ii=VALID_INDEX
    if abs(H_A_list(ii)) > abs(H_B_list(ii))
        v=abs(H_B_list(ii));
    else
        v=abs(H_A_list(ii));
    end
    if v>maxHAB
        maxHAB=v;
    end
end


alpha = 0.5;

global MOD
if strcmp(MOD,'BPSK')
    npairs = 4;
    data_pair= [ 1  1;
            -1  1;
             1 -1;
            -1 -1];
elseif strcmp(MOD,'QPSK')
    npairs = 16;
    data_pair= [ 
             1+j  1+j; 
            1+j  -1+j;
             1+j  1-j;
             1+j  -1-j;
             
             -1+j  1+j;
           -1+j  -1+j;
           -1+j  1-j;
             -1+j  -1-j;
             
             1-j  1+j;
             1-j  -1+j;
             1-j  1-j;
             1-j  -1-j;
             
             -1-j  1+j;
             -1-j  -1+j;
             -1-j  1-j;
             -1-j  -1-j
    ]; 
end
d = zeros(npairs,1);



        for pair_index = 1 : npairs
            d(pair_index) = norm( H_A * data_pair(pair_index, 1)  ...
                + H_B * data_pair(pair_index,2)  ...
                - in )^2;
        end

        
if strcmp(DECODER_TYPE,'hard') == 1
    [value, number] = min(d);
    if strcmp(MOD,'BPSK')
        Sign_data_pair = sign(real(data_pair(number,:)));
        pair = (1+Sign_data_pair)/2;
        xor = bitxor(pair(1),pair(2));
        output_pair = [xor pair];
    end
	
    if strcmp(MOD,'QPSK')
        Sign_data_pair_real = sign(real(data_pair(number,:)));   %qpsk: real and image
        Sign_data_pair_imag = sign(imag(data_pair(number,:)));
        pair_real = (1+Sign_data_pair_real)/2;
        pair_imag = (1+Sign_data_pair_imag)/2;
        xor_real = bitxor(pair_real(1),pair_real(2));
        xor_imag = bitxor(pair_imag(1),pair_imag(2));
        
        output_pair = [xor_real xor_imag pair_real(1) pair_imag(1) pair_real(2) pair_imag(2)];   % re-order at the end
    end
end

if strcmp(DECODER_TYPE, 'soft') == 1
    
  % normal decoder
  if (knownTxRawData(1) == -10) && ((knownTxRawData(2) == -10)) && (knownTxRawData(3) == -10)
  
      if strcmp(MOD,'BPSK')
          %-for user A
          if (d(1)>d(3)) && (d(2)>d(4))
              softA=real(H_A)*real(in+H_B)+imag(H_A)*imag(in+H_B);
          elseif (d(1)>d(3)) && (d(2)<=d(4))
              softA=real(in)*real(H_A-H_B)+imag(in)*imag(H_A-H_B);
          elseif (d(1)<=d(3)) && (d(2)>d(4))
              softA=real(in)*real(H_A+H_B)+imag(in)*imag(H_A+H_B);
          else
              softA=real(H_A)*real(in-H_B)+imag(H_A)*imag(in-H_B);
          end
          NA=maxHA;
          vA = round((((softA/NA)*alpha)/maxHA+0.5)*255);
          if vA<0
              vA=0;
          elseif vA>255
              vA=255;
          end
          
          %-for user B
          if (d(1)>d(2)) && (d(3)>d(4))
              softB=real(H_B)*real(in+H_A)+imag(H_B)*imag(in+H_A);
          elseif (d(1)>d(2)) && (d(3)<=d(4))
              softB=real(in)*real(H_B-H_A)+imag(in)*imag(H_B-H_A);
          elseif (d(1)<=d(2)) && (d(3)>d(4))
              softB=real(in)*real(H_A+H_B)+imag(in)*imag(H_A+H_B);
          else
              softB=real(H_B)*real(in-H_A)+imag(H_B)*imag(in-H_A);
          end
          NB=maxHB;
          vB = round((((softB/NB)*alpha)/maxHB+0.5)*255);
          if vB<0
              vB=0;
          elseif vB>255
              vB=255;
          end
          
          %-for XOR packet
          if (d(1)>d(4)) && (d(2)>d(3))
              softX=real(H_A)*real(in+H_B)+imag(H_A)*imag(in+H_B);
          elseif (d(1)>d(4)) && (d(2)<=d(3))
              softX=real(H_B)*real(in+H_A)+imag(H_B)*imag(in+H_A);
          elseif (d(1)<=d(4)) && (d(2)>d(3))
              softX=real(H_B)*real(H_A-in)+imag(H_B)*imag(H_A-in);
          else
              softX=real(H_A)*real(H_B-in)+imag(H_A)*imag(H_B-in);
          end
          NX=maxHAB;
          vX = round((((softX/NX)*alpha)/maxHAB+0.5)*255);
          if vX<0
              vX=0;
          elseif vX>255
              vX=255;
          end
          
          output_pair=[vX vA vB];
          
      elseif strcmp(MOD,'QPSK')
          %-for user A
          
          
          [valueA_i0, numberA_i0]=min([d(1),d(2),d(3),d(4),d(9),d(10),d(11),d(12)]);%min([d(5),d(6),d(7),d(8),d(13),d(14),d(15),d(16)]);
          [valueA_i1, numberA_i1]=min([d(5),d(6),d(7),d(8),d(13),d(14),d(15),d(16)]);%min([d(1),d(2),d(3),d(4),d(9),d(10),d(11),d(12)]);
          
          [valueA_q0, numberA_q0]=min([d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8)]);%min([d(9),d(10),d(11),d(12),d(13),d(14),d(15),d(16)]);
          [valueA_q1, numberA_q1]=min([d(9),d(10),d(11),d(12),d(13),d(14),d(15),d(16)]);%min([d(1),d(2),d(3),d(4),d(5),d(6),d(7),d(8)]);
          
          softA_i = (valueA_i1-valueA_i0)/8;
          softA_q = (valueA_q1-valueA_q0)/8;
          
          NA=maxHA;
          vA_i = round((((softA_i/NA)*alpha)/maxHA+0.5)*255);
          if vA_i<0
              vA_i=0;
          elseif vA_i>255
              vA_i=255;
          end
          
          vA_q = round((((softA_q/NA)*alpha)/maxHA+0.5)*255);
          if vA_q<0
              vA_q=0;
          elseif vA_q>255
              vA_q=255;
          end
          
          %-for user B
          [valueB_i0, numberB_i0]=min([d(1),d(3),d(5),d(7),d(9),d(11),d(13),d(15)]);%min([d(2),d(4),d(6),d(8),d(10),d(12),d(14),d(16)]);
          [valueB_i1, numberB_i1]=min([d(2),d(4),d(6),d(8),d(10),d(12),d(14),d(16)]);%min([d(1),d(3),d(5),d(7),d(9),d(11),d(13),d(15)]);
          
          [valueB_q0, numberB_q0]=min([d(1),d(2),d(5),d(6),d(9),d(10),d(13),d(14)]);%min([d(3),d(4),d(7),d(8),d(11),d(12),d(15),d(16)]);
          [valueB_q1, numberB_q1]=min([d(3),d(4),d(7),d(8),d(11),d(12),d(15),d(16)]);%min([d(1),d(2),d(5),d(6),d(9),d(10),d(13),d(14)]);
          
          softB_i = (valueB_i1-valueB_i0)/8;
          softB_q = (valueB_q1-valueB_q0)/8;
          
          NB=maxHB;
          vB_i = round((((softB_i/NB)*alpha)/maxHB+0.5)*255);
          if vB_i<0
              vB_i=0;
          elseif vB_i>255
              vB_i=255;
          end
          
          vB_q = round((((softB_q/NB)*alpha)/maxHB+0.5)*255);
          if vB_q<0
              vB_q=0;
          elseif vB_q>255
              vB_q=255;
          end
          
          %-for XOR packet
          [valueX_i0, numberX_i0]=min([d(2),d(4),d(5),d(7),d(10),d(12),d(13),d(15)]);%min([d(1),d(3),d(6),d(8),d(9),d(11),d(14),d(16)]);
          [valueX_i1, numberX_i1]=min([d(1),d(3),d(6),d(8),d(9),d(11),d(14),d(16)]);%min([d(2),d(4),d(5),d(7),d(10),d(12),d(13),d(15)]);
          
          [valueX_q0, numberX_q0]=min([d(3),d(4),d(7),d(8),d(9),d(10),d(13),d(14)]);%min([d(1),d(2),d(5),d(6),d(11),d(12),d(15),d(16)]);
          [valueX_q1, numberX_q1]=min([d(1),d(2),d(5),d(6),d(11),d(12),d(15),d(16)]);%min([d(3),d(4),d(7),d(8),d(9),d(10),d(13),d(14)]);
          
          softX_i = (valueX_i1-valueX_i0)/8;
          softX_q = (valueX_q1-valueX_q0)/8;
          
          NX=maxHAB;
          vX_i = round((((softX_i/NX)*alpha)/maxHAB+0.5)*255);
          if vX_i<0
              vX_i=0;
          elseif vX_i>255
              vX_i=255;
          end
          
          vX_q = round((((softX_q/NX)*alpha)/maxHAB+0.5)*255);
          if vX_q<0
              vX_q=0;
          elseif vX_q>255
              vX_q=255;
          end
          
          output_pair=[vX_i vX_q vA_i vA_q vB_i vB_q];
      end
    
 % SIC decoder
 elseif knownTxRawData(3) ~= -10 % Xor packet
   %fprintf(' known X !!! \n');
   %data_pair= [ 1  1;
   %            -1  1;
   %             1 -1;
   %            -1 -1];
   if knownTxRawData(3) == 0 % pair: [1 4]
       softA=real(in)*real(H_A+H_B)+imag(in)*imag(H_A+H_B);
       softB=real(in)*real(H_A+H_B)+imag(in)*imag(H_A+H_B);
   else % pair: [2 3]
       softA=real(in)*real(H_A-H_B)+imag(in)*imag(H_A-H_B);
       softB=real(in)*real(H_B-H_A)+imag(in)*imag(H_B-H_A);
   end
   
   NA=maxHA;
   vA = round((((softA/NA)*alpha)/maxHA+0.5)*255);
   if vA<0
       vA=0;
   elseif vA>255
       vA=255;
   end
   NB=maxHB;
   vB = round((((softB/NB)*alpha)/maxHB+0.5)*255);
   if vB<0
       vB=0;
   elseif vB>255
       vB=255;
   end
   output_pair=[0 vA vB];

 elseif knownTxRawData(1) ~= -10 % A packet
    if knownTxRawData(1) == 0 % pair: [2 4]
       softB=real(H_B)*real(in+H_A)+imag(H_B)*imag(in+H_A);
       softX=real(H_B)*real(in+H_A)+imag(H_B)*imag(in+H_A);
    else % pair: [1 3]
       softB=real(H_B)*real(in-H_A)+imag(H_B)*imag(in-H_A); 
       softX=real(H_B)*real(H_A-in)+imag(H_B)*imag(H_A-in);
    end
   NB=maxHB;
   vB = round((((softB/NB)*alpha)/maxHB+0.5)*255);
   if vB<0
       vB=0;
   elseif vB>255
       vB=255;
   end
   NX=maxHAB;
   vX = round((((softX/NX)*alpha)/maxHAB+0.5)*255);
   if vX<0
       vX=0;
   elseif vX>255
       vX=255;
   end
    
   output_pair=[vX 0 vB];
    
 elseif knownTxRawData(2) ~= -10 % B packet
     if knownTxRawData(2) == 0 % pair: [3 4]
         softA=real(H_A)*real(in+H_B)+imag(H_A)*imag(in+H_B);
         softX=real(H_A)*real(in+H_B)+imag(H_A)*imag(in+H_B);
     else % pair: [1 2]
         softA=real(H_A)*real(in-H_B)+imag(H_A)*imag(in-H_B);
         softX=real(H_A)*real(H_B-in)+imag(H_A)*imag(H_B-in);
     end
     NA=maxHA;
     vA = round((((softA/NA)*alpha)/maxHA+0.5)*255);
     if vA<0
         vA=0;
     elseif vA>255
         vA=255;
     end
     NX=maxHAB;
     vX = round((((softX/NX)*alpha)/maxHAB+0.5)*255);
     if vX<0
         vX=0;
     elseif vX>255
         vX=255;
     end
    
    output_pair=[vX vA 0];
 end
end

%note: reorder to make it consistent
% output = [A, B, X]
if strcmp(MOD,'BPSK')
    output_pair=[output_pair(2) output_pair(3) output_pair(1)];
elseif strcmp(MOD,'QPSK')
    output_pair=[output_pair(3) output_pair(4) output_pair(5) output_pair(6) output_pair(1) output_pair(2)];
end