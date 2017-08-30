function [mfPrevArray,ParameterVectors,L3Values,ABR,mfPrevIntv...
    demog,ageMthMax,bCulex] = ...
    RunIntvScenarios(Site,AgeLimits,da,MDAInterval,...
    DrugEfficacy,NumYears,MonthlyMDACov,SwitchMonth,...
    ITNCov,IRSCov,VCparams)
   
    IRSCoverages=ones(1,NumYears*12)*ITNCov/100;
    ITNCoverages=ones(1,NumYears*12)*IRSCov/100;
    
    % MultiVecMBR L3Values ParameterVectors ageMthMax bCulex demog kId1 kId2
    % mfPrevArray1  mfPrevArray2 are loaded.
    load(sprintf('ParamVectors%s_v2.mat',Site));
    if ~exist('ABR','var')
        ABR = MultiVecMBR(1,:);
    end
    
    if ~exist('kId1','var')
        kId = 1:length(ParameterVectors(1,:));
    elseif exist('kId1','var') && (kId1(1) > 0)
        kId = unique([kId1;kId2]); clearvars kId1 kId2;
    else
        kId = unique(kId2); clearvars kId1 kId2;
    end
    
    ParameterVectors = ParameterVectors(:,kId);
    L3Values = L3Values(kId,:);
    MultiVecMBR = ABR/12;
    
    mfPrevIntv =...
        Modelling_MDAplusVC(kId,...
        ParameterVectors,L3Values,demog,ageMthMax,da,bCulex,...
        AgeLimits,DrugEfficacy,MDAInterval,NumYears,...
        MultiVecMBR,IRSCoverages,ITNCoverages,MonthlyMDACov,SwitchMonth,VCparams);
   
end
