% Site data: Mf prevalence, blood sample volume, MDA regimen, MDA frequency,...
% MDA coverage, number of years of treatment, vector control, switch year

function [ DokanTofaMf,GwamlarMf,PiapungMf,...
    DokanTofaVol,GwamlarVol,PiapungVol,...
    DokanTofaReg,GwamlarReg,PiapungReg,...
    DokanTofaFreq,GwamlarFreq,PiapungFreq,...
    DokanTofaMDACov,GwamlarMDACov,PiapungMDACov,...
    DokanTofaNumYears,GwamlarNumYears,PiapungNumYears,...
    DokanTofaVC,GwamlarVC,PiapungVC,...
    DokanTofaSwitchYear,GwamlarSwitchYear,PiapungSwitchYear,...
    DokanTofaITNCov,GwamlarITNCov,PiapungITNCov,...
    DokanTofaIRSCov,GwamlarIRSCov,PiapungIRSCov]...
    = PostIntv_data_SSA
%% Post-intervention Mf Data
% enter overall community prevalence in a single line
% 1st column = month of survey; 2nd: Total number of samples; 3rd: Mf +ves;

DokanTofaMf = [(2007-2003)*12+1 151 round(0.013*151);
               (2008-2003)*12+1 158 0;
               (2009-2003)*12+1 223 round(0.004*223)];
GwamlarMf = [(2006-2002)*12+1 240 round(0.121*240);
             (2007-2002)*12+1 128 round(0.016*128);
             (2008-2002)*12+1 100 round(0.05*100);
             (2009-2002)*12+1 143 round(0.049*143)];
PiapungMf = [(2007-2002)*12+1 187 round(0.096*187);
             (2009-2002)*12+1 291 round(0.021*291)];

%% Blood sample volume used to test presence of mf (in uL)
% Standard volumes are 1000 uL, 100 uL, or 20 uL (60 uL uses 100 uL
% correction)
% Model operates based on 1 mL samples, so smaller samples will be
% corrected to reflect this

DokanTofaVol = 60; 
GwamlarVol = 60;
PiapungVol = 60;

%% MDA Regimen 
% if more than one treatment regimen, given frequency for each regimen
% (i.e. [5,2] for ALB followed by DEC+ALB)

% 1: IVM 
% 2: DEC+ALB
% 3: IVM+DEC+ALB
% 4: IVM+ALB
% 5: ALB
% 6: DEC salt

DokanTofaReg = [4]; 
GwamlarReg = [4];
PiapungReg = [4];

%% MDA Frequency
% in months
% if more than one treatment regimen, given frequency for each regimen
% (i.e. [12,6] for annual followed by biannial)

DokanTofaFreq = [12]; 
GwamlarFreq = [12];
PiapungFreq = [12];

%% Annual MDA Coverage
% enter as proportion (i.e. 0.8)

DokanTofaMDACov = [0.83 0.81 0.77 0.9 0.87 0.87 0.87]; % 2003-2009
GwamlarMDACov = [0.92 0.81 0.79 0.75 0.9 0.89 0.91 0.91]; % 2002-2009
PiapungMDACov = [0.85 0.81 0.89 0.96 0.9 0.9 0.9]; % 2003-2009

%% Total number of years of treatment
% (last year - first year) + 1
DokanTofaNumYears = (2009-2003)+1; 
GwamlarNumYears = (2009-2002)+1; 
PiapungNumYears = (2009-2003)+1; 

%% Vector control
% 0: no VC, 1: VC

DokanTofaVC = 0; 
GwamlarVC = 0;
PiapungVC = 0;

%% Year to switch treatment regimens
% enter 0 if a single treatment regimen is followed

DokanTofaSwitchYear = 0; 
GwamlarSwitchYear = 0;
PiapungSwitchYear = 0;

%% Vector control Coverage
% enter as proportion (i.e. 0.8)

DokanTofaITNCov = 0; 
GwamlarITNCov = 0;
PiapungITNCov = 0;

DokanTofaIRSCov = 0; 
GwamlarIRSCov = 0;
PiapungIRSCov = 0;

end
