% min and max values for uniform prior distributions of model parameters 
function [rangeParamVals,maxParamVals,minParamVals] = Range_of_Parameters

% beta;alpha;k0;kLin;k1;r1;sigma1;psi1;psi2s2;mu;gamma;b1;c;ageLev;VoverH;k2;gam2;immC;slopeC;PP;del;alpha2;gamma2
minParamVals=[ 5;0.25;0.000036;0.00000024;3;0.04;1.5;0.1;0.00003;0.008;0.08;0.251;0.015;20;10;3;0.04;0.5;0.01;1;0.001;0.5;0.0125];
maxParamVals=[15;1.50;0.000775;0.28200000;5;0.25;8.5;0.8;0.00364;0.018;0.12;0.485;0.025;30;10;5;0.25;5.5;0.20;9;0.010;9.370;0.2];

rangeParamVals=maxParamVals-minParamVals;
end
