%% add CommonFunctions folder containing basic model functions to path
addpath('../CommonFunctions/');

%% User-defined inputs

% Parameters related to sampling and resampling
NParamVecs = 100000; % default: 200000
SIR_samples = 500; % default: 500

% Site details: Relevant endemic regions, Country (pertains only Sub-Saharan Africa cases), Village
Regions = {'SSA';}; % at least one required
Countries = {'Tanzania'; 'Kenya'; 'Tanzania'}; % required only for SSA sites, names should follow format of Countries_w_Demo.mat
SSAsites = {'Masaika'; 'Kingwede'; 'Kirare'};

% Site data: Mf prevalence, blood sample volume, ABR, vector genera
% arranged in region-specific .m data files created by the user
[ MasaikaMf,KingwedeMf,KirareMf,...
    MasaikaVol,KingwedeVol,KirareVol,...
    MasaikaABR,KingwedeABR,KirareABR,...
    MasaikabCulex,KingwedebCulex,KirarebCulex] = baseline_data_SSA;

% Number of workers for parallel pool settings
NumWorkers = 3; % dependent on system capabilities

%% Import country-level age-demographics of SSA countries if needed

if ~isempty(Countries)
    Import_subSaharanCountries_demoA_demoB;
    load('Countries_w_demo.mat','Country_demo');
end

%% Set up parallel pool
p = gcp('nocreate');
if isempty(p)
    cc=parcluster();
    t=tempname();
    mkdir(t);
    cc.JobStorageLocation=t;
    if exist('parpool') % >= 2013b
        parpool(cc, NumWorkers);
    else % < 2013b
        matlabpool(cc, NumWorkers);
    end
end

%% Run fitting procedure for each site

for iReg = 1:length(Regions)
    Sites  = eval(sprintf('%ssites',Regions{iReg}));
    for iSites = 1%:length(Sites)
        % define data
        MfData  = eval(sprintf('%sMf',Sites{iSites}));
        Vol     = eval(sprintf('%sVol',Sites{iSites}));
        bCulex  = eval(sprintf('%sbCulex',Sites{iSites}));
        ABR     = eval(sprintf('%sABR',Sites{iSites}));
        
        % define demographic parameters
        [demoA1, demoB1] = get_demoA_and_demoB(Countries{iSites},...
            Country_demo,demoA,demoB);
        
        % correct mf data for blood sample volume
        if Vol == 20
            MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.95));
        elseif Vol == 100
            MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.15));
        end
        
        % find maximum age in population
        if isnan(MfData)
            ageMax = 69;
        else
            ageMax = max(MfData(:,4));
        end
        
        % run appropriate fitting procedure
        RunSelectParameterVectors_AgeMf_asInput_wCFA(Sites{iSites},NParamVecs,...
            MfData,ABR,SIR_samples,demoA1,demoB1,...
            ageMax,bCulex);
        
    end
end