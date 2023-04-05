%% Generate SIC phy maps
%  @author: You Lizhao
%  @note: To merge two results to get SIC maps. Here we use exising XOR packet output result.

%InputFiles = {'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
%InputFiles = {'0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB'};
%InputFiles = {'0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB'};
InputFiles = {'0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB','0131_10_74dB_123dB'};
InputAlpha = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};

for fileSICIndex=1:length(InputFiles)
    tag_infile=InputFiles{fileSICIndex};
    alpha = InputAlpha{fileSICIndex};

    ss=strcat('matlab_data/sic_phy_logs/',tag_infile);
    ss=strcat(ss,'_softSIC_alpha_0.05_map.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    
    ss = strcat('matlab_data/rmud_phy/',tag_infile);
    %ss = strcat(ss,'_softNCMA_map.mat');
    ss = strcat(ss,'_maps.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    
    softPhySIC1XorMap = softPhySIC1Map;
    softPhySIC1XorMap(:,X_INDEX) = phyRawMap(:,X_INDEX); 
    
    softPhySIC2XorMap = softPhySIC2Map;
    softPhySIC2XorMap(:,X_INDEX) = phyRawMap(:,X_INDEX); 
    
    ss=strcat('matlab_data/',tag_infile);
    ss=strcat(ss,'_softSIC_map.mat');
    fprintf(' \n Writing file: %s \n', ss);
    save(ss,'alpha','softPhyRawMap','softPhySIC1Map','softPhySIC1XorMap','softPhySIC2Map','softPhySIC2XorMap','multiUser');
end
%}