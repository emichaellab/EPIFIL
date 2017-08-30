function RunSelectParameterVectors_AgeMf_asInput_wCFA(Sites,...
    NParamVecs,MfData,ABR,SIR_samples,demoA,...
    demoB,ageMax,bCulex)

tic; % start timing the run

ageMthMax = ageMax*12; % convert max age in years to months 
demog=[1, 1, demoA, demoB, 1]; % set up demographic parameters
da=1; % integration step-size of 1 month
toleranceX=0.000001; % Used in endemic equilibrium of state variables

% Set random number stream
RandStream.setGlobalStream...
    (RandStream('mt19937ar','seed',rand*sum(10000*clock)));

% Initialize Output Arrays
ParameterVectors = [];
L3Values1 = [];
mfPrevArray1 = [];
cfaPrevArray1 = [];
ABR1 = [];

% Loop to generate an initial sample of parameter vectors, run the model
% to equilibrium, and resample best fits until the desired number of best
% fits defined by SIR_samples is achieved
NumParam = 0;
while ( NumParam < SIR_samples )
    
    % generate an initial sample of parameter vectors
    ParamVectors = ParameterVectors_LHSDesign_based(NParamVecs);
    ParamVectors(15,:) = (ABR/12)./ParamVectors(1,:); % VoverH = MBR/beta 
    
    % run the model to equilibrium using the sampled parameter vectors
    [L3Values,mfPrevArray,cfaPrevArray] = ...
        Calculate_EndemicEquil_L3_Mf_CFA(NParamVecs,...
        ageMthMax,bCulex,demog,ParamVectors,da,toleranceX);
    
    % Resample the parameter vectors using 1) SIR algorithm and 2) pass/fail criteria
    [kId1, ~] = get_modelFits_Mf_AgeData(MfData,...
        mfPrevArray,demog,da,ageMthMax,SIR_samples);
    kId = unique(kId1); 
    
    % Store parameters, equilibrium L3, and an equilibrium mf values
    % corresponding to the resampled parameter vectors 
    ParameterVectors(:,NumParam+1:NumParam+length(kId)) = ...
        ParamVectors(:,kId);
    
    L3Values1(NumParam+1:NumParam+length(kId),:) = ...
        L3Values(kId,:);
    
    mfPrevArray1(:,NumParam+1:NumParam+length(kId)) = ...
        mfPrevArray(:,kId);
    
    cfaPrevArray1(:,NumParam+1:NumParam+length(kId)) = ...
        cfaPrevArray(:,kId);
    
    ABR1(NumParam+1:NumParam+length(kId),:) = ...
        12*ParamVectors(15,kId)'.*ParamVectors(1,kId)';
    
    % Update number of parameter vectors selected 
    NumParam = NumParam + length(kId);
    fprintf(1,'Num selected kIds = %d, Total Num Params = %d\n',...
        length(kId),NumParam);
end

L3Values = L3Values1;
mfPrevArray = mfPrevArray1;
cfaPrevArray = cfaPrevArray1;
ABR = ABR1;
toc; % stop timing the run

% save sampled best fitting outputs and basic descriptive parameters of
% the site
save(sprintf('ParamVectors%s.mat',Sites),'bCulex',...
    'demog','ABR','ParameterVectors','L3Values','ageMthMax',...
    'mfPrevArray','cfaPrevArray');

end
% #########################################################################
function [kId1, kId2] = get_modelFits_Mf_AgeData(MfData,...
    mfPrevArray,demog,da,ageMthMax,SIR_samples)

% Initialize the index and likelihood arrays
kId1 = []; % for indices of SIR resampled fits
kId2 = []; % for indices of pass/fail criteria fits
LikArray1 = ones(length(mfPrevArray(1,:)),1); % for likelihoods of SIR resampled fits
LikArray2 = zeros(length(mfPrevArray(1,:)),1); % for likelihoods of pass/fail criteria fits

% calculate 95% CI upper and lower bounds of each mf data point
MfBounds = get_the95LU_bounds_agedata(MfData);

% calculate likelihood of each sampled parameter vector using SIR and
% pass/fail methods
parfor i=1:length(mfPrevArray(1,:))
    LikArray1(i)= calculateLikelihoods(MfData,100*mfPrevArray(:,i),demog,da,ageMthMax);
    LikArray2(i)= likelihood_using_ranks_by_passFail_Mf_Only(mfPrevArray(:,i),MfBounds);
end

% resample by SIR algorithm
Indx1 = ReSampleMfPrevalence(LikArray1,SIR_samples);
Indx = find( Indx1 > 0);
if ~isempty(Indx)
    kId1 = [kId1; Indx1(Indx)];
end

% resample by pass/fail criteria
Indx2 = ReSampleMfPrevalence_passFail_Mf_only(LikArray2,SIR_samples);
Indx = find( Indx2 > 0);
if ~isempty(Indx)
    kId2 = [kId2; Indx2(Indx)];
end

end
% #########################################################################
