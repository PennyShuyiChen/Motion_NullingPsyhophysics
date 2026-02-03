%% Run ANOVA analyses demo 
% Modified -Penny-Shuyi Chen 
% Feb 03, 2026
%% 1. Organize data (one value per trial)
%  -y
 trialData.YpeakCorr = PkMaxsY;
 trialData.YpeakLag = PkLagsY;
 trialData.YFWHM = FWHM_Y;
 trialData.YFWHMgauss = FWHMgauss_Y;
 for k=1:size(condSFStrVec,2)
 trialData.SF{k} = char(condSFStrVec(k));
 end
 trialData.eccentricity = str2double(condStrVec);
 for k=1:size(subVec,2)
 trialData.subject{k} = char(subVec(k));
 end
 %% -x 
 trialData.XpeakCorr = PkMaxsX;
 trialData.XpeakLag = PkLagsX;
 trialData.XFWHM = FWHM_X;
 trialData.XFWHMgauss = FWHMgauss_X;
 for k=1:size(condSFStrVec,2)
 trialData.SF{k} = char(condSFStrVec(k));
 end
 trialData.eccentricity = str2double(condStrVec);
 for k=1:size(subVec,2)
 trialData.subject{k} = char(subVec(k));
 end
 
 %% -phase 
 trialData.PhpeakCorr = PkMaxsPh;
 trialData.PhpeakLag = PkLagsPh;
 trialData.PhFWHM = FWHM_Ph;
 trialData.PhFWHMgauss = FWHMgauss_Ph;
 for k=1:size(condSFStrVec,2)
 trialData.SF{k} = char(condSFStrVec(k));
 end
 trialData.eccentricity = str2double(condStrVec);
 for k=1:size(subVec,2)
 trialData.subject{k} = char(subVec(k));
 end

%% 2. Define factors and DVs
% factorNames = {'trialType', 'duration'};
% dvNames = {'peakCorr', 'peakLag'};
% dvLabels = {'Peak Correlation', 'Peak Lag (s)'};
factorNames = {'SF', 'eccentricity'};
dvNames = {'YpeakCorr', 'YpeakLag','YFWHM','YFWHMgauss'};
dvLabels = {'Peak Correlation', 'Peak Lag (s)','FWHM','Guassian FWHM'};

%% 3. Run ANOVAs
clc;
results = runCCGAnova(trialData, factorNames, dvNames, dvLabels);

%% 4. Plot
plotCCGAnovaResults(results, trialData);
