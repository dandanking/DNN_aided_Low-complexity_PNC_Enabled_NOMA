
function [output_pair,debug] = SIC_Subcarrier_Demod(in,H_A,H_B,H_A_list,H_B_list,knownTxRawData)
global DECODER_TYPE
global alpha;

%{
VALID_INDEX = [7:32 34:59];
minHA = min(abs(H_A_list(VALID_INDEX)));
minHB = min(abs(H_B_list(VALID_INDEX)));
minHAB = 1000;
maxHAB = 0;
for ii=VALID_INDEX
    if abs(H_A_list(ii)) > abs(H_B_list(ii))
        v=abs(H_B_list(ii));
    else
        v=abs(H_A_list(ii));
    end
    if v<minHAB
        minHAB=v;
    end
    
    v = abs(H_B_list(ii))+abs(H_A_list(ii));
    if v>maxHAB
        maxHAB=v;
    end
end
%}

%fprintf(' In SIC Decoder \n');

A_INDEX=1; B_INDEX=2; INVALID=-10;
output_pair=[INVALID,INVALID,INVALID];
debug=[0,0];

if strcmp(DECODER_TYPE,'hard') == 1
    if (knownTxRawData(A_INDEX) == -10) && ((knownTxRawData(B_INDEX) == -10))
        
        symA=in/H_A;
        dataA=sign(real(symA));
        bitsA=(1+dataA)/2;
        
        symB=in/H_B;
        dataB=sign(real(symB));
        bitsB=(1+dataB)/2;
        
        output_pair=[bitsA,bitsB,INVALID];
        debug=[symA,symB];
        
    elseif knownTxRawData(A_INDEX) ~= -10 % A packet
        if knownTxRawData(A_INDEX) == 0
            inCancel = in - H_A*(-1);
        else
            inCancel = in - H_A*(1);
        end
        symB = inCancel/H_B;
        dataB = sign(real(symB));
        bitsB = (1+dataB)/2;
        
        output_pair=[0 bitsB INVALID];
        debug=[0,symB];
        
    elseif knownTxRawData(B_INDEX) ~= -10 % B packet
        if knownTxRawData(B_INDEX) == 0
            inCancel = in - H_B*(-1);
        else
            inCancel = in - H_B*(1);
        end
        symA = inCancel/H_A;
        dataA = sign(real(symA));
        bitsA = (1+dataA)/2;
        
        output_pair=[bitsA 0 INVALID];
        debug=[symA,0];
    end
end

maxHAB = abs(H_A)+abs(H_B);
if strcmp(DECODER_TYPE,'soft') == 1
    if (knownTxRawData(A_INDEX) == -10) && ((knownTxRawData(B_INDEX) == -10))
        
        %symA=in*abs(H_A)/(H_A*maxHAB);
        %dataA=round((real(symA)*alpha+1.5)*255);
        %symB=in*abs(H_B)/(H_B*maxHAB);
        %dataB=round((real(symB)*alpha+1.5)*255);
        symA = in/H_A; symB = in/H_B;
        dataA=clamp(real(in/H_A),alpha,8);
        dataB=clamp(real(in/H_B),alpha,8);
        
        output_pair=[dataA,dataB,INVALID];
        debug=[symA,symB];
        
    elseif knownTxRawData(A_INDEX) ~= -10 % A packet
        if knownTxRawData(A_INDEX) == 0
            inCancel = in - H_A*(-1);
        else
            inCancel = in - H_A*(1);
        end
        %symB = inCancel*abs(H_B)/(H_B*maxHAB);
        %dataB = round((real(symB)*alpha+0.5)*255);
        symB = inCancel/H_B;
        dataB = clamp(real(inCancel/H_B),alpha,8);
        
        output_pair=[0 dataB INVALID];
        debug=[0,symB];
        
    elseif knownTxRawData(B_INDEX) ~= -10 % B packet
        if knownTxRawData(B_INDEX) == 0
            inCancel = in - H_B*(-1);
        else
            inCancel = in - H_B*(1);
        end
        %symA = inCancel*abs(H_A)/(H_A*maxHAB);
        %dataA = round((real(symA)*alpha+0.5)*255);
        symA = inCancel/H_A;
        dataA = clamp(real(inCancel/H_A),alpha,8);
        
        output_pair=[dataA 0 INVALID];
        debug=[symA,0];
    end
end


%{
alpha = 0.5;

d = zeros(4,1);

data_pair= [ 1  1;
            -1  1;
             1 -1;
            -1 -1];

        for pair_index = 1 : 4;
            d(pair_index) = norm( H_A * data_pair(pair_index, 1)  ...
                + H_B * data_pair(pair_index,2)  ...
                - in )^2;
        end

if strcmp(DECODER_TYPE, 'soft') == 1
    
  if (knownTxRawData(A_INDEX) == -10) && ((knownTxRawData(B_INDEX) == -10))
  
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

  end

end

%}
%note: reorder to make it consistent
%output_pair=[output_pair(2) output_pair(3) output_pair(1)];
