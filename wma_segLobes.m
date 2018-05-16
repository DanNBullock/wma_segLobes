function [classification] =wma_segLobes(feORwbfg, fsDir)
%
% [classification] =wma_segLobes(wbfg, fsDir)
% This function automatedly segments the intralobal fibers of the brain

% Inputs:
% -wbfg: a whole brain fiber group structure
% -fsDir: path to THIS SUBJECT'S freesurfer directory

% Outputs:
% classification:  a classification structure with .name and .indexes
% fields
%
% (C) Daniel Bullock, 2018, Indiana University

%% parameter note & initialization


[wbfg, fe] = bsc_LoadAndParseFiberStructure(feORwbfg)

labelNifti= wma_getAsegFile(fsDir , '2009');
LeftStreams=midpoints(:,1)<0;

classification=[];

classification.index=zeros(length(wbfg.fibers),1);
classification.names={'Left_Frontal','Right_Frontal','Left_Parietal','Right_Parietal','Left_Temporal','Right_Temporal','Left_Occipital','Right_Occipital','Left_Cerebellar','Right_Cerebellar'};

cerebellumROI=[8 7; 46 47];






for leftright= [1,2]
    
    
    sidenum=10000+leftright*1000;  

       
    
    FrontalROI=bsc_roiFromFSnums(fsDir,[124 148 118 165 101 154 105 115 154 155 115 170 129 146 153 ...
        164 106 116 108 131 171 112 150 104 169 114 113]+sidenum ,0);
    [Frontal, FrontalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {FrontalROI}], { 'both_endpoints' }, 'Frontal');
    
    TemporalROI=bsc_roiFromFSnums(fsDir,[144 134 138 137 173 174 135 175 121 151 123 162]+sidenum ,0);
    [Temporal, TemporalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {TemporalROI}], { 'both_endpoints' }, 'Temporal');
    
    occipitalROI=bsc_roiFromFSnums(fsDir,[120 119 111 158 166 143 145 159 152 122 162 161 121]+sidenum,0);
    [Occipital, OccipitalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {occipitalROI}], { 'both_endpoints' }, 'Occipital');
    
    ParietalROI=bsc_roiFromFSnums(fsDir,[157 127 168 136 126 125 156 128 141 172 147 109]+sidenum,0);
    [Parietal, ParietalIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {ParietalROI}], { 'both_endpoints' }, 'Parietal');
    
    CerebellarROI=bsc_roiFromFSnums(fsDir,cerebellumROI(leftright,:),0);
    [Cerebellar, CerebellarIND]=wma_SegmentFascicleFromConnectome(wbfg, [  {CerebellarROI}], { 'both_endpoints' }, 'Parietal');
    
    
    
    
    if leftright==1
        fprintf('\n Left segmentation complete')
        classification.index(FrontalIND & LeftStreams)=find( strcmp(classification.names,'Left Frontal'));
        classification.index(TemporalIND & LeftStreams)=find( strcmp(classification.names,'Left Temporal'));
        classification.index(OccipitalIND & LeftStreams)=find( strcmp(classification.names,'Left Occipital'));
        classification.index(ParietalIND & LeftStreams)=find( strcmp(classification.names,'Left Temporal'));
        classification.index(CerebellarIND & LeftStreams)=find( strcmp(classification.names,'Left Cerebellar'));
    else
        fprintf('\n Right segmentation complete')
        classification.index(FrontalIND & ~LeftStreams)=find( strcmp(classification.names,'Right Frontal'));
        classification.index(TemporalIND & ~LeftStreams)=find( strcmp(classification.names,'Right Temporal'));
        classification.index(OccipitalIND & ~LeftStreams)=find( strcmp(classification.names,'Right Occipital'));
        classification.index(ParietalIND & ~LeftStreams)=find( strcmp(classification.names,'Right Temporal'));
        classification.index(CerebellarIND & ~LeftStreams)=find( strcmp(classification.names,'Right Cerebellar'));
    end
    
end
end