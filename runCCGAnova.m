%% CCG ANOVA FUNCTIONS
% Functions for running trial-level and repeated measures ANOVAs
% on cross-correlation (CCG) data from psychophysics experiments.
%
% These functions can be adapted for different experiments with varying:
%   - Number of subjects
%   - Number and types of conditions (factors)
%   - Different dependent variables
%
% USAGE:
%   1. Organize your trial data with CCG metrics and condition labels
%   2. Call runCCGAnova() to compute both trial-level and repeated measures ANOVAs
%   3. Call plotCCGAnovaResults() to visualize the results
%
% Author: Sufia Ahmad
% Date: March 2024
% Modified Penny-Shuyi Chen Feb 03, 2026


%%  FUNCTION 1: runCCGAnova
%Runs both trial-level and repeated measures ANOVAs on CCG data
function results = runCCGAnova(trialData, factorNames, dvNames, dvLabels)

% INPUTS:
%   trialData   - Structure with fields:
%                   .dv1, .dv2, ... (dependent variables, e.g., peakCorr, peakLag)
%                   .factor1, .factor2, ... (condition labels for each trial)
%                   .subject (subject ID for each trial - REQUIRED for repeated measures)
%                 All fields should be vectors of the same length (one element per trial)
%
%   factorNames - Cell array of factor field names in trialData
%                 e.g., {'trialType', 'duration', 'contrast'}
%
%   dvNames     - Cell array of dependent variable field names
%                 e.g., {'peakCorr', 'peakLag', 'width'}
%
%   dvLabels    - Cell array of labels for plotting (same order as dvNames)
%                 e.g., {'Peak Correlation', 'Peak Lag (ms)', 'CCG Width (ms)'}
%
% OUTPUT:
%   results     - Structure containing:
%                   .trialLevel  - Trial-level ANOVA results (each trial = observation)
%                   .subjAvg     - Subject-averaged repeated measures ANOVA results
%                   .factorNames - Factor names used
%                   .dvNames     - DV names used
%                   .dvLabels    - DV labels
%                   .nTrials     - Number of trials
%                   .nSubjects   - Number of subjects
%                   .subjects    - List of unique subjects
%
% EXAMPLE:
%   % Prepare data
%   trialData.peakCorr = [0.25, 0.28, 0.22, ...];  % one value per trial
%   trialData.peakLag = [45, 52, 38, ...];
%   trialData.trialType = {'continuous', 'continuous', 'delayed', ...};
%   trialData.duration = [1, 1, 2, ...];
%   trialData.contrast = [1, 2, 1, ...];
%   trialData.subject = {'S1', 'S1', 'S2', ...};
%
%   % Run ANOVAs
%   results = runCCGAnova(trialData, {'trialType', 'duration', 'contrast'}, ...
%                         {'peakCorr', 'peakLag'}, {'Peak Correlation', 'Peak Lag (ms)'});

%% Validate inputs
if nargin < 4
    error('Usage: runCCGAnova(trialData, factorNames, dvNames, dvLabels)');
end

if ~isfield(trialData, 'subject')
    error('trialData must have a "subject" field for repeated measures ANOVA');
end

% Check all required fields exist
for i = 1:length(factorNames)
    if ~isfield(trialData, factorNames{i})
        error('trialData is missing factor field: %s', factorNames{i});
    end
end
for i = 1:length(dvNames)
    if ~isfield(trialData, dvNames{i})
        error('trialData is missing DV field: %s', dvNames{i});
    end
end

%% Initialize results structure
results = struct();
results.factorNames = factorNames;
results.dvNames = dvNames;
results.dvLabels = dvLabels;

%% Get basic info
nTrials = length(trialData.(dvNames{1}));
subjects = unique(trialData.subject);
nSubjects = length(subjects);

results.nTrials = nTrials;
results.nSubjects = nSubjects;
results.subjects = subjects;

%% Prepare factor arrays
nFactors = length(factorNames);
factors = cell(1, nFactors);
for fi = 1:nFactors
    factors{fi} = trialData.(factorNames{fi});
    % Convert to categorical if cell array
    if iscell(factors{fi})
        factors{fi} = categorical(factors{fi});
    else
        factors{fi} = categorical(factors{fi});
    end
end
subjectCat = categorical(trialData.subject);

%% Trial level ANOVA
fprintf('\n--- TRIAL-LEVEL ANOVA (each trial as independent observation) ---\n');
fprintf('WARNING: inflates N and may have inflated Type I error rate.\n\n');

results.trialLevel = struct();

for dvi = 1:length(dvNames)
    dvName = dvNames{dvi};
    dvLabel = dvLabels{dvi};
    dv = trialData.(dvName);
    
    % Remove NaN values
    validIdx = ~isnan(dv);
    dvValid = dv(validIdx);
    factorsValid = cell(1, nFactors);
    for fi = 1:nFactors
        factorsValid{fi} = factors{fi}(validIdx);
    end
    
    fprintf('%s:\n', dvLabel);
    
    try
        % Run N-way ANOVA
        [p, tbl, stats] = anovan(dvValid, factorsValid, ...
            'model', 'full', ...
            'varnames', factorNames, ...
            'display', 'off');
        
        % Calculate effect sizes (partial eta-squared)
        SS_error = tbl{end, 2};
        effectSizes = zeros(length(p), 1);
        for i = 1:length(p)
            SS_effect = tbl{i+1, 2};
            effectSizes(i) = SS_effect / (SS_effect + SS_error);
        end
        
        % Store results
        results.trialLevel.(dvName).p = p;
        results.trialLevel.(dvName).table = tbl;
        results.trialLevel.(dvName).stats = stats;
        results.trialLevel.(dvName).effectSizes = effectSizes;
        results.trialLevel.(dvName).label = dvLabel;
        results.trialLevel.(dvName).nObs = length(dvValid);
        
        % Print results
        printAnovaTable(tbl, p, effectSizes, factorNames);
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results.trialLevel.(dvName).error = ME.message;
    end
end

%% Subject averaged repeated measures ANOVA
fprintf('\n--- SUBJECT-AVERAGED ANOVA (repeated measures, N = %d subjects) ---\n', nSubjects);
fprintf('approach for repeated measures designs.\n\n');

results.subjAvg = struct();

% Get unique levels of each factor
factorLevels = cell(1, nFactors);
for fi = 1:nFactors
    factorLevels{fi} = unique(factors{fi});
end

for dvi = 1:length(dvNames)
    dvName = dvNames{dvi};
    dvLabel = dvLabels{dvi};
    dv = trialData.(dvName);
    
    fprintf('%s:\n', dvLabel);
    
    % Compute subject means for each condition combination
    [subjMeans, subjFactors, subjIDs] = computeSubjectMeans(dv, factors, ...
        subjectCat, subjects, factorLevels, factorNames);
    
    if isempty(subjMeans)
        fprintf('  ERROR: Could not compute subject means\n');
        results.subjAvg.(dvName).error = 'Could not compute subject means';
        continue;
    end
    
    % Remove NaN
    validIdx = ~isnan(subjMeans);
    subjMeans = subjMeans(validIdx);
    for fi = 1:nFactors
        subjFactors{fi} = subjFactors{fi}(validIdx);
    end
    subjIDs = subjIDs(validIdx);
    
    try
        % Build model matrix for anovan
        % Include all main effects and interactions except those involving subject
        nFactorsWithSubj = nFactors + 1;  % factors + subject
        
        % Run ANOVA with subject as random factor
        allFactors = [subjFactors, {categorical(subjIDs)}];
        allVarnames = [factorNames, {'Subject'}];
        
        [p, tbl, stats] = anovan(subjMeans, allFactors, ...
            'model', 'interaction', ...
            'random', nFactorsWithSubj, ...  % Subject is random
            'varnames', allVarnames, ...
            'display', 'off');
        
        % Calculate effect sizes
        SS_error = tbl{end, 2};
        effectSizes = zeros(length(p), 1);
        for i = 1:min(length(p), size(tbl,1)-2)
            SS_effect = tbl{i+1, 2};
            if ~isnan(SS_effect) && ~isnan(SS_error) && SS_error > 0
                effectSizes(i) = SS_effect / (SS_effect + SS_error);
            end
        end
        
        % Store results
        results.subjAvg.(dvName).p = p;
        results.subjAvg.(dvName).table = tbl;
        results.subjAvg.(dvName).stats = stats;
        results.subjAvg.(dvName).effectSizes = effectSizes;
        results.subjAvg.(dvName).label = dvLabel;
        results.subjAvg.(dvName).nSubjects = nSubjects;
        results.subjAvg.(dvName).nObs = length(subjMeans);
        results.subjAvg.(dvName).subjMeans = subjMeans;
        results.subjAvg.(dvName).subjFactors = subjFactors;
        results.subjAvg.(dvName).subjIDs = subjIDs;
        
        % Print results (skip subject effect)
        printAnovaTableRM(tbl, p, effectSizes, allVarnames);
        
    catch ME
        fprintf('  ERROR: %s\n', ME.message);
        results.subjAvg.(dvName).error = ME.message;
    end
end



end

%% FUNCTION 2: plotCCGAnovaResults
%  Plots ANOVA results (main effects and interactions)

function plotCCGAnovaResults(results, trialData, varargin)
% INPUTS:
%   results     - Output from runCCGAnova()
%   trialData   - Same trialData structure used in runCCGAnova()
%   
%   Optional Name-Value pairs:
%     'plotType'    - 'both' (default), 'trialLevel', or 'subjAvg'
%     'colors'      - Structure with color definitions (optional)
%
% EXAMPLE:
%   plotCCGAnovaResults(results, trialData);
%   plotCCGAnovaResults(results, trialData, 'plotType', 'subjAvg');

%% Parse inputs
p = inputParser;
addRequired(p, 'results');
addRequired(p, 'trialData');
addParameter(p, 'plotType', 'both', @ischar);
addParameter(p, 'colors', [], @isstruct);
parse(p, results, trialData, varargin{:});

plotType = p.Results.plotType;

%% Get info from results
factorNames = results.factorNames;
dvNames = results.dvNames;
dvLabels = results.dvLabels;
nFactors = length(factorNames);
nDVs = length(dvNames);
nSubjects = results.nSubjects;

%% Compute data for plotting
% Get unique levels of each factor
factorLevels = cell(1, nFactors);
factorLabels = cell(1, nFactors);
for fi = 1:nFactors
    vals = trialData.(factorNames{fi});
    if iscell(vals)
        factorLevels{fi} = unique(vals);
    else
        factorLevels{fi} = unique(vals);
    end
    factorLabels{fi} = cellstr(string(factorLevels{fi}));
end

%% Generate default colors
nColors = max(cellfun(@length, factorLevels));
defaultColors = lines(nColors);

%% Trial level plots
if strcmp(plotType, 'both') || strcmp(plotType, 'trialLevel')
    
    figure('Position', [50, 50, 400*nFactors, 300*nDVs], ...
        'Name', 'ANOVA Main Effects - Trial Level');
    
    for dvi = 1:nDVs
        dvName = dvNames{dvi};
        dv = trialData.(dvName);
        
        for fi = 1:nFactors
            subplot(nDVs, nFactors, (dvi-1)*nFactors + fi);
            hold on;
            
            factor = trialData.(factorNames{fi});
            levels = factorLevels{fi};
            nLevels = length(levels);
            
            % Compute means and SEMs for each level
            means = nan(1, nLevels);
            sems = nan(1, nLevels);
            
            for li = 1:nLevels
                if iscell(factor)
                    idx = strcmp(factor, levels{li});
                else
                    idx = factor == levels(li);
                end
                vals = dv(idx);
                vals = vals(~isnan(vals));
                means(li) = mean(vals);
                sems(li) = std(vals) / sqrt(length(vals));
            end
            
            % Plot bars
            b = bar(1:nLevels, means, 0.6);
            b.FaceColor = 'flat';
            for li = 1:nLevels
                b.CData(li,:) = defaultColors(li,:);
            end
            
            % Error bars
            errorbar(1:nLevels, means, sems, 'k', 'LineStyle', 'none', ...
                'LineWidth', 1.5, 'CapSize', 8);
            
            % Add significance
            if isfield(results.trialLevel, dvName) && ~isfield(results.trialLevel.(dvName), 'error')
                pVal = results.trialLevel.(dvName).p(fi);
                addSignificance(means, sems, pVal);
            end
            
            % Labels
            set(gca, 'XTick', 1:nLevels, 'XTickLabel', factorLabels{fi}, 'FontSize', 10);
            if fi == 1
                ylabel(dvLabels{dvi}, 'FontSize', 11);
            end
            if dvi == 1
                title(factorNames{fi}, 'FontSize', 12, 'FontWeight', 'bold');
            end
            set(gca, 'Box', 'off', 'TickDir', 'out');
        end
    end
    
    sgtitle('Main Effects - Trial Level ANOVA', 'FontSize', 14, 'FontWeight', 'bold');
end

%% Subject average plots 
if strcmp(plotType, 'both') || strcmp(plotType, 'subjAvg')
    
    figure('Position', [100, 100, 400*nFactors, 300*nDVs], ...
        'Name', 'ANOVA Main Effects - Subject Averaged (Repeated Measures)');
    
    % First compute subject-level means for each condition
    subjects = results.subjects;
    
    for dvi = 1:nDVs
        dvName = dvNames{dvi};
        dv = trialData.(dvName);
        
        for fi = 1:nFactors
            subplot(nDVs, nFactors, (dvi-1)*nFactors + fi);
            hold on;
            
            factor = trialData.(factorNames{fi});
            levels = factorLevels{fi};
            nLevels = length(levels);
            
            % Compute subject means for each level, then average across subjects
            means = nan(1, nLevels);
            sems = nan(1, nLevels);
            
            for li = 1:nLevels
                subjMeansForLevel = nan(nSubjects, 1);
                
                for si = 1:nSubjects
                    % Find trials for this subject and level
                    if iscell(factor)
                        levelIdx = strcmp(factor, levels{li});
                    else
                        levelIdx = factor == levels(li);
                    end
                    
                    if iscell(trialData.subject)
                        subjIdx = strcmp(trialData.subject, subjects{si});
                    else
                        subjIdx = trialData.subject == subjects(si);
                    end
                    
                    idx = levelIdx & subjIdx;
                    vals = dv(idx);
                    vals = vals(~isnan(vals));
                    
                    if ~isempty(vals)
                        subjMeansForLevel(si) = mean(vals);
                    end
                end
                
                % Average across subjects
                validSubj = ~isnan(subjMeansForLevel);
                means(li) = mean(subjMeansForLevel(validSubj));
                sems(li) = std(subjMeansForLevel(validSubj)) / sqrt(sum(validSubj));
            end
            
            % Plot bars
            b = bar(1:nLevels, means, 0.6);
            b.FaceColor = 'flat';
            for li = 1:nLevels
                b.CData(li,:) = defaultColors(li,:);
            end
            
            % Error bars
            errorbar(1:nLevels, means, sems, 'k', 'LineStyle', 'none', ...
                'LineWidth', 1.5, 'CapSize', 8);
            
            % Add significance
            if isfield(results.subjAvg, dvName) && ~isfield(results.subjAvg.(dvName), 'error')
                pVal = results.subjAvg.(dvName).p(fi);
                addSignificance(means, sems, pVal);
            end
            
            % Labels
            set(gca, 'XTick', 1:nLevels, 'XTickLabel', factorLabels{fi}, 'FontSize', 10);
            if fi == 1
                ylabel(dvLabels{dvi}, 'FontSize', 11);
            end
            if dvi == 1
                title(factorNames{fi}, 'FontSize', 12, 'FontWeight', 'bold');
            end
            set(gca, 'Box', 'off', 'TickDir', 'out');
        end
    end
    
    sgtitle(sprintf('Main Effects - Repeated Measures ANOVA (N = %d subjects)', nSubjects), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

%% Interaction plot (for 2+ factors) 
if nFactors >= 2 && (strcmp(plotType, 'both') || strcmp(plotType, 'subjAvg'))
    
    figure('Position', [150, 150, 400*nDVs, 350], ...
        'Name', 'Interaction Effects - Subject Averaged');
    
    % Plot Factor1 x Factor2 interaction
    factor1 = trialData.(factorNames{1});
    factor2 = trialData.(factorNames{2});
    levels1 = factorLevels{1};
    levels2 = factorLevels{2};
    nLevels1 = length(levels1);
    nLevels2 = length(levels2);
    
    colors1 = lines(nLevels1);
    
    for dvi = 1:nDVs
        subplot(1, nDVs, dvi);
        hold on;
        
        dvName = dvNames{dvi};
        dv = trialData.(dvName);
        
        % Compute means for each combination
        plotHandles = [];
        legendLabels = {};
        
        for l1i = 1:nLevels1
            means = nan(1, nLevels2);
            sems = nan(1, nLevels2);
            
            for l2i = 1:nLevels2
                subjMeans = nan(nSubjects, 1);
                
                for si = 1:nSubjects
                    % Find trials for this combination
                    if iscell(factor1)
                        idx1 = strcmp(factor1, levels1{l1i});
                    else
                        idx1 = factor1 == levels1(l1i);
                    end
                    if iscell(factor2)
                        idx2 = strcmp(factor2, levels2{l2i});
                    else
                        idx2 = factor2 == levels2(l2i);
                    end
                    if iscell(trialData.subject)
                        idxS = strcmp(trialData.subject, subjects{si});
                    else
                        idxS = trialData.subject == subjects(si);
                    end
                    
                    idx = idx1 & idx2 & idxS;
                    vals = dv(idx);
                    vals = vals(~isnan(vals));
                    
                    if ~isempty(vals)
                        subjMeans(si) = mean(vals);
                    end
                end
                
                validSubj = ~isnan(subjMeans);
                means(l2i) = mean(subjMeans(validSubj));
                sems(l2i) = std(subjMeans(validSubj)) / sqrt(sum(validSubj));
            end
            
            % Plot line with error bars
            h = errorbar((1:nLevels2) + (l1i-1)*0.1 - 0.05, means, sems, '-o', ...
                'Color', colors1(l1i,:), 'LineWidth', 2, ...
                'MarkerSize', 8, 'MarkerFaceColor', colors1(l1i,:), 'CapSize', 6);
            plotHandles(end+1) = h;
            legendLabels{end+1} = string(levels1(l1i));
        end
        
        % Add interaction significance
        if isfield(results.subjAvg, dvName) && ~isfield(results.subjAvg.(dvName), 'error')
            % Find interaction p-value (usually index nFactors+1 for 2-way)
            if length(results.subjAvg.(dvName).p) > nFactors
                pVal = results.subjAvg.(dvName).p(nFactors + 1);
                yLims = ylim;
                if pVal < 0.001
                    sigStr = 'Interaction: ***';
                elseif pVal < 0.01
                    sigStr = 'Interaction: **';
                elseif pVal < 0.05
                    sigStr = 'Interaction: *';
                else
                    sigStr = 'Interaction: n.s.';
                end
                text(mean(1:nLevels2), yLims(2)*0.95, sigStr, ...
                    'HorizontalAlignment', 'center', 'FontSize', 11, 'FontWeight', 'bold');
            end
        end
        
        set(gca, 'XTick', 1:nLevels2, 'XTickLabel', factorLabels{2}, 'FontSize', 10);
        xlabel(factorNames{2}, 'FontSize', 11);
        ylabel(dvLabels{dvi}, 'FontSize', 11);
        
        if dvi == 1
            legend(plotHandles, legendLabels, 'Location', 'best', 'EdgeColor', 'none');
        end
        
        title(sprintf('%s × %s', factorNames{1}, factorNames{2}), ...
            'FontSize', 12, 'FontWeight', 'bold');
        set(gca, 'Box', 'off', 'TickDir', 'out');
        xlim([0.5, nLevels2 + 0.5]);
    end
    
    sgtitle(sprintf('Interaction Effects - Repeated Measures (N = %d)', nSubjects), ...
        'FontSize', 14, 'FontWeight', 'bold');
end

end

%% HELPER FUNCTIONS

function [subjMeans, subjFactors, subjIDs] = computeSubjectMeans(dv, factors, ...
    subjectCat, subjects, factorLevels, factorNames)
% Compute subject means for each condition combination

nFactors = length(factors);
nSubjects = length(subjects);

% Generate all combinations of factor levels
nLevelsPerFactor = cellfun(@length, factorLevels);
nCombinations = prod(nLevelsPerFactor);

% Initialize outputs
subjMeans = [];
subjFactors = cell(1, nFactors);
for fi = 1:nFactors
    subjFactors{fi} = {};
end
subjIDs = {};

% Loop through all subjects and all condition combinations
for si = 1:nSubjects
    subj = subjects(si);
    subjIdx = subjectCat == subj;
    
    % Generate all combinations using ndgrid
    grids = cell(1, nFactors);
    [grids{:}] = ndgrid(factorLevels{:});
    
    for ci = 1:nCombinations
        % Get the factor levels for this combination
        condIdx = true(size(dv));
        condIdx = condIdx & subjIdx;
        
        for fi = 1:nFactors
            levelVal = grids{fi}(ci);
            condIdx = condIdx & (factors{fi} == levelVal);
        end
        
        % Compute mean for this subject in this condition
        vals = dv(condIdx);
        vals = vals(~isnan(vals));
        
        if ~isempty(vals)
            subjMeans(end+1) = mean(vals);
            for fi = 1:nFactors
                subjFactors{fi}{end+1} = char(grids{fi}(ci));
            end
            subjIDs{end+1} = char(subj);
        end
    end
end

end

function printAnovaTable(tbl, p, effectSizes, factorNames)
% Print formatted ANOVA table

% Generate effect names (main effects and interactions)
effectNames = generateEffectNames(factorNames);

fprintf('  %-30s %8s %8s %10s %6s\n', 'Effect', 'F', 'p', 'η²p', 'Sig');
fprintf('  %-30s %8s %8s %10s %6s\n', '------------------------------', '--------', '--------', '----------', '------');

for i = 1:min(length(p), length(effectNames))
    F = tbl{i+1, 6};
    pVal = p(i);
    eta = effectSizes(i);
    
    [pStr, sigStr] = formatPValue(pVal);
    
    fprintf('  %-30s %8.2f %8s %10.3f %6s\n', effectNames{i}, F, pStr, eta, sigStr);
end
fprintf('\n');

end

function printAnovaTableRM(tbl, p, effectSizes, varnames)
% Print formatted ANOVA table for repeated measures (skip Subject effect)

fprintf('  %-30s %8s %8s %10s %6s\n', 'Effect', 'F', 'p', 'η²p', 'Sig');
fprintf('  %-30s %8s %8s %10s %6s\n', '------------------------------', '--------', '--------', '----------', '------');

for i = 1:length(p)
    % Skip if this is the Subject effect
    effectName = tbl{i+1, 1};
    if contains(effectName, 'Subject')
        continue;
    end
    
    F = tbl{i+1, 6};
    pVal = p(i);
    if i <= length(effectSizes)
        eta = effectSizes(i);
    else
        eta = NaN;
    end
    
    [pStr, sigStr] = formatPValue(pVal);
    
    fprintf('  %-30s %8.2f %8s %10.3f %6s\n', effectName, F, pStr, eta, sigStr);
end
fprintf('\n');

end

function effectNames = generateEffectNames(factorNames)
% Generate effect names including main effects and interactions

nFactors = length(factorNames);
effectNames = factorNames;  % Start with main effects

% Add 2-way interactions
if nFactors >= 2
    for i = 1:nFactors
        for j = i+1:nFactors
            effectNames{end+1} = sprintf('%s × %s', factorNames{i}, factorNames{j});
        end
    end
end

% Add 3-way interaction
if nFactors >= 3
    effectNames{end+1} = sprintf('%s × %s × %s', factorNames{1}, factorNames{2}, factorNames{3});
end

end

function [pStr, sigStr] = formatPValue(pVal)
% Format p-value and significance string

if pVal < 0.001
    pStr = '<.001';
    sigStr = '***';
elseif pVal < 0.01
    pStr = sprintf('%.3f', pVal);
    sigStr = '**';
elseif pVal < 0.05
    pStr = sprintf('%.3f', pVal);
    sigStr = '*';
else
    pStr = sprintf('%.3f', pVal);
    sigStr = '';
end

end

function addSignificance(means, sems, pVal)
% Add significance marker to bar plot

yMax = max(means + sems);
yMin = min(means - sems);
yRange = yMax - yMin;
sigY = yMax + yRange * 0.1;

nBars = length(means);
plot([1, nBars], [sigY, sigY], 'k-', 'LineWidth', 1.5);

if pVal < 0.001
    sigStr = '***';
elseif pVal < 0.01
    sigStr = '**';
elseif pVal < 0.05
    sigStr = '*';
else
    sigStr = 'n.s.';
end

text(mean(1:nBars), sigY + yRange * 0.05, sigStr, ...
    'HorizontalAlignment', 'center', 'FontSize', 14, 'FontWeight', 'bold');

end
