%% add CommonFunctions folder containing basic model functions and Baseline folder containing baseline fits to path
addpath('../CommonFunctions/');
addpath('../BaselineFitting/');

%% User-defined inputs

% integration time step
da = 1; % 1 month

% Site details: Relevant endemic regions, Country (pertains only Sub-Saharan Africa cases), Village
Regions = {'SSA';}; % at least one required
Countries = {'Nigeria';'Nigeria';'Nigeria'}; % required only for SSA sites, names should follow format of Countries_w_Demo.mat
SSAsites = {'DokanTofa'; 'Gwamlar'; 'Piapung';};

% Site data: Mf prevalence, blood sample volume, MDA regimen, MDA frequency,...
% MDA coverage, number of years of treatment, vector control, switch year
% arranged in region-specific .m data files created by the user
[ DokanTofaMf,GwamlarMf,PiapungMf,...
    DokanTofaVol,GwamlarVol,PiapungVol,...
    DokanTofaReg,GwamlarReg,PiapungReg,...
    DokanTofaFreq,GwamlarFreq,PiapungFreq,...
    DokanTofaMDACov,GwamlarMDACov,PiapungMDACov,...
    DokanTofaNumYears,GwamlarNumYears,PiapungNumYears,...
    DokanTofaVC,GwamlarVC,PiapungVC,...
    DokanTofaSwitchYear,GwamlarSwitchYear,PiapungSwitchYear,...
    DokanTofaITNCov,GwamlarITNCov,PiapungITNCov,...
    DokanTofaIRSCov,GwamlarIRSCov,PiapungIRSCov]...
    = PostIntv_data_SSA;

% Number of workers for parallel pool settings
NumWorkers = 3; % dependent on system capabilities

% Specify if CFA should be calculated
CFAFlag = 1;

%% Intervention specifications

% Mass drug administration parameters
AgeLimits=[5 100]; % treatment for > 5 yrs old
RegimenEfficacy0 = [0.1 0.99 9; % IVM worm kill rate, mf kill rate, sterilization period (months)
    0.55 0.95 6; % DEC+ALB
    0.55 1 120; % IVM+DEC+ALB (permanent sterilization - 10 yrs)
    0.35 0.99 9; % IVM+ALB
    0.55 0.95 6; % ALB
    0.59 0.86 10]; % DEC salt

% Vector control parameters
AnnualDecrease = 0;
IRSParams = [87.5, 93, 277, 6, 80]; % Vector control by indoor residual spraying
ITNParams = [20, 90, 97, 26, 12*3, 80]; % Vector control by long lasting insecticide nets
VCparams = SSA_IRS_ITN_Parameters(IRSParams,ITNParams,AnnualDecrease);

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
        %% set up inputs
        
        % define data
        MfData  = eval(sprintf('%sMf',Sites{iSites}));
        Vol     = eval(sprintf('%sVol',Sites{iSites}));
        MDAReg = eval(sprintf('%sReg',Sites{iSites}));
        MDAFreq = eval(sprintf('%sFreq',Sites{iSites}));
        MDACov  = eval(sprintf('%sMDACov',Sites{iSites}));
        VC  = eval(sprintf('%sVC',Sites{iSites}));
        IRSCov  = eval(sprintf('%sIRSCov',Sites{iSites}));
        ITNCov  = eval(sprintf('%sITNCov',Sites{iSites}));
        NumYears  = eval(sprintf('%sNumYears',Sites{iSites}));
        SwitchYear  = eval(sprintf('%sSwitchYear',Sites{iSites}));
        
        % convert data into format compatible with model
        
        % correct mf data for blood sample volume
        if ~isnan(MfData(1,1)) && (Vol == 20)
            MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.95));
        elseif ~isnan(MfData(1,1)) && (Vol == 100 || Vol == 60)
            MfData(:,3)=floor(min(MfData(:,2),MfData(:,3)*1.15));
        end
        
        % Month to switch drug treatment regimen
        SwitchMonth = zeros(1,length(SwitchYear));
        for h = 1:length(SwitchYear)
            if SwitchYear(h) == 0
                SwitchMonth(h) = 0;
            else
                SwitchMonth(h) = SwitchYear(h)*12;
            end
        end
        
        % MDA Frequency
        MDAInterval = zeros(1,NumYears*12);
        if length(MDAFreq)==1
            MDAInterval(1,1:end) = MDAFreq(1);
        else
            MDAInterval(1,1:SwitchMonth) = MDAFreq(1);
            MDAInterval(1,SwitchMonth+1:end) = MDAFreq(2);
        end
        
        % MDA Drug Regimen Efficacy
        RegimenEfficacy = zeros(length(MDAReg),length(RegimenEfficacy0(1,:)));
        for j = 1:length(MDAReg)
            RegimenEfficacy(j,:) = RegimenEfficacy0(MDAReg(j),:);
        end
        
        % MDA Coverage
        MonthlyMDACov = zeros(NumYears*12,1);
        for i = 1:length(MDACov)
            MonthlyMDACov((i-1)*MDAFreq(1)+1:i*MDAFreq(1)) = MDACov(i);
        end
        
        % VC coverages
        IRSCoverages=ones(1,NumYears*12)*ITNCov;
        ITNCoverages=ones(1,NumYears*12)*IRSCov;
        
        %% load baseline fits
        % ABR, L3Values, ParameterVectors, ageMthMax, bCulex, demog, mfPrevArray are loaded.
        load(sprintf('ParamVectors%s.mat',Sites{iSites}));
        MultiVecMBR = ABR/12;
        kId = 1:length(ABR);
                      
        %% model interventions
        [mfPrevIntv,MBRIntv,L3Intv] =...
            Modelling_MDAplusVC(kId,...
            ParameterVectors,L3Values,demog,ageMthMax,da,bCulex,...
            AgeLimits,RegimenEfficacy,MDAInterval,NumYears,...
            MultiVecMBR,IRSCoverages,ITNCoverages,MonthlyMDACov,SwitchMonth,VCparams);
        
        save(sprintf('Intv%s.mat',char(Sites{iSites})),...
                'mfPrevIntv','MBRIntv','L3Intv','RegimenEfficacy','MonthlyMDACov','IRSCoverages',...
                'ITNCoverages','MfData','SwitchMonth');
    end
end
