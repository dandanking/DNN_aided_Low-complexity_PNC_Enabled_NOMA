
function v = clamp(in,alpha,nbits)
base = 2^(nbits-1);
v = round((in*alpha+1)*base);
if v < 0
    v = 0;
end
if v > 2^nbits-1
    v = nbits-1;
end
end