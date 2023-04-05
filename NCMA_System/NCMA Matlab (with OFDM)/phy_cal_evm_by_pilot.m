function [ evm ] = phy_cal_evm_by_pilot( measure_energy,nsymbol)
measure_energy=measure_energy.';

evm=[];

T=2*nsymbol;  % 2 pilots for both A and B in one symbol
Aideal=1;     % for BPSK
Cideal=1;     % pilots use 1

for num=1:length(measure_energy)/T
 
  Vmeas=measure_energy( ( (num-1)*T+1): num*T);
  Pv=0;              %average power in one pkt
  for i=1:T
      Pv=Pv+abs(Vmeas(i))^2;
  end

  Ameas=sqrt(T/Pv);


  temp=0;
 
  for j=1:length(T)
      temp=temp+(real(Vmeas(j))*Ameas-real(Cideal)*Aideal)^2+(imag(Vmeas(j))*Ameas-imag(Cideal)*Aideal)^2;
  end

  evm=[evm,sqrt(temp/T)];
 
 end

end

