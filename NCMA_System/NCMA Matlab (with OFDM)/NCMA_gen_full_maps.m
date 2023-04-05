%% Generate full phy maps
%  @author: You Lizhao
%  @note: To merge PMUD, SIC2, and XOR maps
%  @note: full maps means PMUD + SIC2 + XOR maps

%{
files = {'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
for fileIndex=1:length(files)
    files = {'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
    tag=files{fileIndex};

    ss=strcat('matlab_data/sic_tmp/',tag);
    ss=strcat(ss,'_softSIC_alpha_0.05_map.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    map1 = softPhySIC2Map;
    
    ss = strcat('matlab_data/',tag);
    ss = strcat(ss,'_maps.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    softPhyFullMap =  phyRawMap;

    for ii=1:length(softPhyFullMap)
        if  map1(ii,A_INDEX) == 1
            softPhyFullMap(ii,A_INDEX) = 1;
        end
        if  map1(ii,B_INDEX) == 1
            softPhyFullMap(ii,B_INDEX) = 1;
        end
    end
    
    ss=strcat('matlab_data/',tag);
    ss=strcat(ss,'_softFull_map.mat');
    save(ss,'softPhyFullMap','multiUser');
end
%}

A_INDEX = 1; B_INDEX = 2;
%InputFiles = {'0128_011_75dB','0126_09_8dB','0125_08_85dB','0128_01_9dB','0126_01_95dB','0126_04_10dB','0126_06_105dB'};
%InputFiles = {'0130_03_95dB_75dB','0130_04_95dB_85dB','0130_05_95dB_95dB','0130_06_95dB_105dB','0130_08_95dB_135dB'};
%InputFiles = {'0128_11_75dB_75dB','0128_013_75dB_95dB','0130_01_75dB_105dB','0128_14_75dB_135dB'};
InputFiles = {'0131_01_01_209dB_123dB','0131_02_01_90dB_123dB','0131_07_01_70dB_74dB','0131_08_01_70dB_90dB','0131_10_74dB_123dB'};
files={'0131_07_01_70dB_74dB'};
InputAlpha = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
for fileIndex=1:length(InputFiles)
    tag   = InputFiles{fileIndex};
    alpha = InputAlpha{fileIndex};

    ss=strcat('matlab_data/sic_phy/',tag);
    ss=strcat(ss,'_softSIC_map.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    map1 = softPhySIC2Map;
    
    ss = strcat('matlab_data/rmud_phy/',tag);
    %ss = strcat(ss,'_softNCMA_map.mat');
    ss = strcat(ss,'_maps.mat');
    fprintf(' \n Opening file: %s \n', ss);
    load(ss);
    
    softPhyFullMap =  phyRawMap;

    for ii=1:length(softPhyFullMap)
        if  map1(ii,A_INDEX) == 1
            softPhyFullMap(ii,A_INDEX) = 1;
        end
        if  map1(ii,B_INDEX) == 1
            softPhyFullMap(ii,B_INDEX) = 1;
        end
    end
    
    ss=strcat('matlab_data/',tag);
    ss=strcat(ss,'_softFull_map.mat');
    save(ss,'softPhyFullMap','multiUser');
    
end
%}
