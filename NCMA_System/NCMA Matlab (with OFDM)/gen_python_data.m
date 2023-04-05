

clear all;

files={'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
%files={'0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB'};
%files={'0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB'};
%files={'0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB'};

for fileIndex=1:length(files)

tag=files{fileIndex};
ss=strcat('matlab_data/rmud_phy/',tag);
ss_in=strcat(ss,'_maps.mat');
fprintf(' Opening file: %s \n', ss_in);
load(ss_in);

ss=strcat('python_data/rmud_phy/',tag);
ss_out1=strcat(ss,'_raw_map.dat');
ss_out2=strcat(ss,'_bridge_map.dat');
ss_out3=strcat(ss,'_twophase_map.dat');

dlmwrite(ss_out1,phyRawMap);
dlmwrite(ss_out2,phyBridgeMap);
dlmwrite(ss_out3,phyTwoPhaseMap);


end
