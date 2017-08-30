% Sample baseline data for sites Masaika,TZ, Kingwede,KEN, and Kirare,TZ 

function [ MasaikaMf,KingwedeMf,KirareMf,...
    MasaikaVol,KingwedeVol,KirareVol,...
    MasaikaABR,KingwedeABR,KirareABR,...
    MasaikabCulex,KingwedebCulex,KirarebCulex] = baseline_data_SSA

%% Baseline Mf Data
% Best data to provide would be age-stratified prevalence
% Variable names = {SiteName}Mf

% 1st column = mid-age of group; 2nd: Total number of samples; 3rd: Mf +ves;
% 4th: upper age of group

% If age-stratified data is not available, enter overall community
% prevalence in a single line following the column definitions above

MasaikaMf = [...
    2.5,   84,  1,  4
    7.0,  103,  8,  9
    12.0, 129, 23, 14
    17.0,  80, 20, 19
    24.5, 142, 47, 29
    34.5, 122, 34, 39
    44.5,  77, 29, 49
    54.5,  41, 19, 59
    64.5,  70, 30, 69];

KingwedeMf = [...
    2.5,  141, 0,  4
    7.0,  139, 0,  9
    12.0, 145, 1, 14
    17.0,  78, 2, 19
    24.5, 125, 6, 29
    34.5,  86, 8, 39
    44.5,  37, 1, 49
    54.5,  41, 1, 59
    64.5,  33, 3, 69];

KirareMf = [...
    5	324	23	10
    15	251	58	20
    25	148	42	30
    35	137	37	40
    45	89	21	50
    55	55	15	60
    65	67	22	70];

%% Blood sample volume used to test presence of mf (in uL)
% Standard volumes are 1000 uL, 100 uL, or 20 uL
% Model operates based on 1 mL samples, so smaller samples will be
% corrected to reflect this

MasaikaVol = 1000;
KingwedeVol = 1000;
KirareVol = 1000;

%% Baseline ABR data
% If not available for a particular site, enter value as NaN

MasaikaABR = 6184;
KingwedeABR = 1548;
KirareABR = 2091;

%% Mosquito species flag based on dominant genera (0: Anopheles, 1: Culex)

MasaikabCulex = 0;
KingwedebCulex = 1;
KirarebCulex = 0;

end
