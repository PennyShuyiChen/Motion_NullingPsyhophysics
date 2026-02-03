% doit - modified for the hinge approximation for V1 cells 
% PSC ver. 2026/02/02 
    % NOTE check the layout in ccgSF for both dependencies

%%
close all; clear all;clc;
%% trkDat class
% Functions
% trkDat('filename')
t = trkDat('cfcvp.mat');   % files you can test with 'zebv.mat', 'cfv.mat'                         
t.makeccgs;

%% ccg class - diff SFs *** Finalized for both SF and ecc dependence - ccgSF ***
%rd = x.plotraw_sf('c');

 x = ccgSF('xc.mat');
 br = x.bootanal_sf('c'); % for individual subject plots 
 %br = x.bootanal_sf('s'); % for popuation plots 
 x.br = br;
 
 %acc = x.Acc_sf;
 
 %% Get the trial data structure ready for repeated-measure ANOVA (per task)
 % compute dependent property Lags (which is in secs) 
 load('cfcLags.mat');% get 'theLags'
 load('cfcEcc.mat')
 load('cfcSF.mat')
 load('cfcSub.mat')
 zeroInd = ceil(length(theLags)./2);
 pLags = theLags(zeroInd:end);
 
 % trial data
 xCorY = x.xc.xCorY(:,zeroInd:end); % y-vertical direction tracking
 xCorX = x.xc.xCorX(:,zeroInd:end); % x-horizontal direction tracking
 xCorPh = -x.xc.xCorDP(:,zeroInd:end);% phase tracking 
 
 %% get the parameters for rmANOVA (eg. pkCor, lag, FWHM)
 PkMaxsY = nan(1,1200); PkLagsY = nan(1,1200); PkIdxY = nan(1,1200);
 PkMaxsX = nan(1,1200); PkLagsX = nan(1,1200); PkIdxX = nan(1,1200);
 PkMaxsPh = nan(1,1200); PkLagsPh = nan(1,1200); PkIdxPh = nan(1,1200);
 gauss = @(b,x) b(1)*exp(-(x-b(2)).^2/(2*b(3)^2));
 
 %%
 for k = 1:size(xCorY,1)
     Y = xCorY(k,:);
     [PkMaxsY(k), PkIdxY(k)] = max(Y);
     PkLagsY(k) = pLags(PkIdxY(k)); 
     hMaxY(k) = PkMaxsY(k)./2;
     i1 = find(Y(1:PkIdxY(k)) <= hMaxY(k), 1, 'last');%left crossing 
     i2 = find(Y(PkIdxY(k):end) <= hMaxY(k), 1, 'first') + PkIdxY(k) - 1;
     if isempty(i1)||isempty(i2)
         FWHM_Y(k) = nan;
     else
         x1 = interp1(Y(i1:i1+1), pLags(i1:i1+1), hMaxY(k));      
         x2 = interp1(Y(i2-1:i2), pLags(i2-1:i2), hMaxY(k));%right corssing 
         FWHM_Y(k) = x2-x1;
     end
     % gauss fit ver of FWHM 
     mu0 = pLags(PkIdxY(k));
     sigma0 = (max(pLags) - min(pLags)) / 6;%rough guess 
     opts = optimset('Display','off');
     b = lsqcurvefit(gauss, [PkMaxsY(k) mu0 sigma0], pLags, Y, [], [], opts);
     A     = b(1);
     mu    = b(2);
     sigma = abs(b(3));
     FWHMgauss_Y(k) = 2*sqrt(2*log(2)) * sigma;
        
 end
 
 
 %%
 for k = 1:size(xCorX,1)
     Y = xCorX(k,:);
     [PkMaxsX(k), PkIdxX(k)] = max(xCorX(k,:));
     PkLagsX(k) = pLags(PkIdxX(k));
     
     hMaxX(k) = PkMaxsX(k)./2;
     i1 = find(Y(1:PkIdxX(k)) <= hMaxX(k), 1, 'last');%left crossing 
     i2 = find(Y(PkIdxX(k):end) <= hMaxX(k), 1, 'first') + PkIdxX(k) - 1;
     if isempty(i1)||isempty(i2)
         FWHM_X(k) = nan;
     else
         x1 = interp1(Y(i1:i1+1), pLags(i1:i1+1), hMaxX(k));      
         x2 = interp1(Y(i2-1:i2), pLags(i2-1:i2), hMaxX(k));%right corssing 
         FWHM_X(k) = x2-x1;
     end
     % gauss fit ver of FWHM 
     mu0 = pLags(PkIdxX(k));
     sigma0 = (max(pLags) - min(pLags)) / 6;%rough guess 
     opts = optimset('Display','off');
     b = lsqcurvefit(gauss, [PkMaxsX(k) mu0 sigma0], pLags, Y, [], [], opts);
     A     = b(1);
     mu    = b(2);
     sigma = abs(b(3));
     FWHMgauss_X(k) = 2*sqrt(2*log(2)) * sigma;
 end 
 
 
 %%
 for k = 1:size(xCorPh,1)
     Y = xCorPh(k,:);
     [PkMaxsPh(k), PkIdxPh(k)] = max(xCorPh(k,:));
     PkLagsPh(k) = pLags(PkIdxPh(k));
     
     hMaxPh(k) = PkMaxsPh(k)./2;
     i1 = find(Y(1:PkIdxPh(k)) <= hMaxPh(k), 1, 'last');%left crossing 
     i2 = find(Y(PkIdxPh(k):end) <= hMaxPh(k), 1, 'first') + PkIdxPh(k) - 1;
     if isempty(i1)||isempty(i2)
         FWHM_Ph(k) = nan;
     else
         x1 = interp1(Y(i1:i1+1), pLags(i1:i1+1), hMaxPh(k));      
         x2 = interp1(Y(i2-1:i2), pLags(i2-1:i2), hMaxPh(k));%right corssing 
         FWHM_Ph(k) = x2-x1;
     end
     % gauss fit ver of FWHM 
     mu0 = pLags(PkIdxPh(k));
     sigma0 = (max(pLags) - min(pLags)) / 6;%rough guess 
     opts = optimset('Display','off');
     b = lsqcurvefit(gauss, [PkMaxsPh(k) mu0 sigma0], pLags, Y, [], [], opts);
     A     = b(1);
     mu    = b(2);
     sigma = abs(b(3));
     FWHMgauss_Ph(k) = 2*sqrt(2*log(2)) * sigma;
         
 end 
 
 
 
 

%% ccg class - diff eccs *** Finalized for both SF and ecc dependence - ccgSF ***
x = ccgSF('xc.mat');

brE = x.bootanal_ecc('c');
%brE = x.bootanal_ecc('c');
x.br = brE;

%% ccg class - analysis *** This is the Finalized ver. for ecc dependece only***
%  x = ccgAnalysis('xc.mat');
% % x.plotraw_an('c');
% % x.plotmeans_an('c');
%  br = x.bootanal_an('c');
%  x.br = br;
%  
%   %acceleration 
%  %acc = x.Acc;
 %% Function Fits - Gaussian + DoG  
 %[gaussFit,gf,check,amp] = GaussFit(br);
 [dog1Fit,dgf,dcheck,damp] = dog1fit(br);
 %% Paramter vs. ecc 
 [gaussLMfit, gaussT] = paramLMfit(gaussFit,check);
 %[dog1LMfit, dog1T, rawdata] = paramLMfit(dog1Fit,dcheck);
 
