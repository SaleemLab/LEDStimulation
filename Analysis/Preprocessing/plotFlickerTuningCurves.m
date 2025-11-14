function plotFlickerTuningCurves(spikeData, varInfo, options)
if ~exist('options','var'), options=struct; end
if ~isfield(options,'subplotDimNames'), options.subplotDimNames = {}; end
if ~isfield(options,'overlayDimNames'), options.overlayDimNames = {}; end
if ~isfield(options,'baseline'), options.baseline = 'none'; end
if ~isfield(options,'plotStyles'), options.plotStyles = {'-', '--', ':', '-.'}; end
if ~isfield(options,'plotIndividualPoints'), options.plotIndividualPoints = false; end
if ~isfield(options,'subplotTitles'), options.subplotTitles = {}; end
if ~isfield(options,'baselineCollapseDims'), options.baselineCollapseDims = {}; end
if ~isfield(options,'valueLabels'), options.valueLabels = struct; end

if ~isfield(options,'xAxisDimName')
    error('The required option "xAxisDimName" was not provided. Please specify which variable to plot on the x-axis in the options struct.');
end

allDimNames = {varInfo.name};
[~, xAxisDim] = ismember(options.xAxisDimName, allDimNames);
if xAxisDim == 0, error('xAxisDimName "%s" not found in varInfo.', options.xAxisDimName); end

[~, subplotDims] = ismember(options.subplotDimNames, allDimNames);
if any(subplotDims == 0)
    invalidIdx = find(subplotDims == 0, 1);
    error('subplotDimNames entry "%s" not found in varInfo.', options.subplotDimNames{invalidIdx});
end

[~, overlayDims] = ismember(options.overlayDimNames, allDimNames);
if any(overlayDims == 0)
    invalidIdx = find(overlayDims == 0, 1);
    error('overlayDimNames entry "%s" not found in varInfo.', options.overlayDimNames{invalidIdx});
end

[~, baselineCollapseDimIdx] = ismember(options.baselineCollapseDims, allDimNames);
if any(baselineCollapseDimIdx == 0)
    invalidIdx = find(baselineCollapseDimIdx == 0, 1);
    error('baselineCollapseDims entry "%s" not found in varInfo.', options.baselineCollapseDims{invalidIdx});
end


otherDims = setdiff(1:length(allDimNames), [xAxisDim, subplotDims, overlayDims]);

nSubplotVars = length(subplotDims);
subplotGridSize = [varInfo(subplotDims).nVals];
if isempty(subplotGridSize), subplotGridSize = 1; end

% Handle any number of subplot dimensions by creating an optimal 2D grid
numTotalSubplots = prod(subplotGridSize);
if numTotalSubplots > 0
    % Find the integer factors of numTotalSubplots that are closest to the square root
    nRows = floor(sqrt(numTotalSubplots));
    while mod(numTotalSubplots, nRows) ~= 0
        nRows = nRows - 1;
    end
    nCols = numTotalSubplots / nRows;
else
    nRows = 1;
    nCols = 1;
end


figure('WindowStyle', 'normal', 'Position', [100, 100, 400*nCols, 350*nRows]);

subplotIdxCombs = allcombs(subplotGridSize);
allAxes = []; 

for iSubplot = 1:size(subplotIdxCombs, 1)
    ax = subplot(nRows, nCols, iSubplot);
    allAxes(end+1) = ax;
    hold on;
    
    titleStr = {};
    
    nOverlayVars = length(overlayDims);
    overlayGridSize = [varInfo(overlayDims).nVals];
    if isempty(overlayGridSize), overlayGridSize = 1; end
    overlayIdxCombs = allcombs(overlayGridSize);
    
    for iOverlay = 1:size(overlayIdxCombs, 1)
        
        idx = repmat({':'}, 1, length(allDimNames));
        
        for k = 1:nSubplotVars
            idx{subplotDims(k)} = subplotIdxCombs(iSubplot, k);
            currentVarName = varInfo(subplotDims(k)).name;
            currentVal = varInfo(subplotDims(k)).uniqueVals(subplotIdxCombs(iSubplot, k));
            
            % Check for a custom label for the current value
            displayVal = num2str(currentVal); % Default to the numerical value
            if isfield(options.valueLabels, currentVarName)
                labelMap = options.valueLabels.(currentVarName);
                if isa(labelMap, 'containers.Map') && isKey(labelMap, currentVal)
                    displayVal = labelMap(currentVal); % Use custom string label
                end
            end
            titleStr{k} = sprintf('%s: %s', currentVarName, displayVal);
        end
        
        for k = 1:nOverlayVars
            idx{overlayDims(k)} = overlayIdxCombs(iOverlay, k);
        end
        
        dataSlice = squeeze(spikeData(idx{:}));
        
        xVals = varInfo(xAxisDim).uniqueVals;
        nXPoints = length(xVals);
        
        plotBaseline = false;
        baselineIdx = [];
        
        if isfield(options, 'baseline')
            if ischar(options.baseline)
                if strcmpi(options.baseline, 'last')
                    baselineIdx = nXPoints;
                    plotBaseline = true;
                elseif strcmpi(options.baseline, 'first')
                    baselineIdx = 1;
                    plotBaseline = true;
                elseif ~strcmpi(options.baseline, 'none')
                    warning('Invalid baseline string option. Must be ''first'', ''last'', or ''none''. Ignoring.');
                end
            elseif isnumeric(options.baseline) && isscalar(options.baseline) ...
                    && options.baseline >= 1 && options.baseline <= nXPoints
                baselineIdx = round(options.baseline);
                plotBaseline = true;
            elseif ~isequal(options.baseline, 'none')
                 warning('Invalid baseline option. Must be ''first'', ''last'', an integer index, or ''none''. Ignoring.');
            end
        end

        if plotBaseline
            if ~isempty(options.baselineCollapseDims)
                baselineFetchIdx = repmat({':'}, 1, length(allDimNames));
                baselineFetchIdx{xAxisDim} = baselineIdx;
                
                % Fix overlay dimensions
                for k = 1:nOverlayVars
                    baselineFetchIdx{overlayDims(k)} = overlayIdxCombs(iOverlay, k);
                end
                
                % Fix subplot dimensions that are NOT being collapsed
                for k = 1:nSubplotVars
                    if ~ismember(subplotDims(k), baselineCollapseDimIdx)
                        baselineFetchIdx{subplotDims(k)} = subplotIdxCombs(iSubplot, k);
                    end
                end
                
                baselineDataCells = spikeData(baselineFetchIdx{:});
                baselineData = cat(2, baselineDataCells{:});

            else
                baselineData = dataSlice{baselineIdx};
            end
            
            tuningIdx = true(1, nXPoints);
            tuningIdx(baselineIdx) = false;
            tuningData = dataSlice(tuningIdx);
            
            xTickLabels = xVals(tuningIdx);
            xPlotPoints = 1:length(xTickLabels);
            
            mean_baseline = mean(mean(baselineData, 1, 'omitnan'), 2, 'omitnan');
            sem_baseline = sem(mean(baselineData, 1, 'omitnan'), 2);
            shadedErrorBar(xPlotPoints, repelem(mean_baseline, 1, length(xPlotPoints)), ...
                           repelem(sem_baseline, 1, length(xPlotPoints)), ...
                           'lineProps', {options.plotStyles{iOverlay}, 'LineStyle', ':'}, 'patchSaturation', 0.1);
        else
            tuningData = dataSlice;
            xTickLabels = xVals;
            xPlotPoints = 1:length(xTickLabels);
        end
        
        mean_rates = cellfun(@(x) mean(mean(x, 1, 'omitnan'), 2, 'omitnan'), tuningData);
        
        if options.plotIndividualPoints
            % Plot the mean line without error bars
            h = plot(xPlotPoints, mean_rates, options.plotStyles{iOverlay});
            
            % Overlay individual trial points
            lineColor = h.Color;
            for iPoint = 1:numel(tuningData)
                trialData = tuningData{iPoint};
                if ~isempty(trialData)
                    trialMeans = mean(trialData, 1, 'omitnan');
                    x_jitter = (rand(size(trialMeans)) - 0.5) * 0.25;
                    x_pos = xPlotPoints(iPoint) + x_jitter;
                    plot(x_pos, trialMeans, 'o', 'Color', lineColor, 'MarkerFaceColor', lineColor, 'MarkerSize', 4, 'LineStyle', 'none');
                end
            end
        else
            % Default behavior: plot with error bars
            sem_rates = cellfun(@(x) sem(mean(x, 1, 'omitnan'), 2), tuningData);
            errorbar(xPlotPoints, mean_rates, sem_rates, options.plotStyles{iOverlay});
        end
    end
    
    hold off;
    
    if ~isempty(options.subplotTitles)
        if numel(options.subplotTitles) >= iSubplot
            finalTitleStr = options.subplotTitles{iSubplot};
        else
            warning('Not enough subplotTitles provided. Reverting to default for subplot %d.', iSubplot);
            finalTitleStr = strjoin(titleStr, ', ');
        end
    else
        finalTitleStr = strjoin(titleStr, ', ');
    end
    title(finalTitleStr, 'Interpreter', 'none');

    xlabel(varInfo(xAxisDim).name, 'Interpreter', 'none');
    ylabel('Mean Spike Count');
    ax = gca;
    ax.XTick = xPlotPoints;
    ax.XTickLabel = xTickLabels;
    xlim([min(xPlotPoints)-0.5, max(xPlotPoints)+0.5]);
    defaultAxesProperties(gca,true)
end

linkaxes(allAxes, 'y');

end

function s = sem(data, dim)
    if nargin < 2, dim = find(size(data) ~= 1, 1); end
    if isempty(dim), dim = 1; end
    n = sum(~isnan(data), dim);
    s = std(data, 0, dim, 'omitnan') ./ sqrt(n);
end

function combs = allcombs(v)
    if isempty(v)
        combs = [];
        return;
    end
    n = length(v);
    gridArgs = arrayfun(@(x) 1:x, v, 'UniformOutput', false);
    V = cell(1, n);
    [V{:}] = ndgrid(gridArgs{:});
    combs = reshape(cat(n+1, V{:}), [], n);
end

function defaultAxesProperties(ax, offsetFlag)
set(ax, 'TickDir','out', 'TickLength', [0.01 0.001]...
    ,  'color', 'none', 'box','off', 'XColor', 'k', 'YColor', 'k',...
    'FontName', 'Calibri', 'LineWidth', 0.5)
if nargin>1
    if offsetFlag
        offsetAxes(ax);
    end
end
end

function offsetAxes(ax)
    % This is a helper function for defaultAxesProperties.
    % The main visual effect of offsetting axes is achieved by setting
    % 'box' to 'off' and 'TickDir' to 'out' in defaultAxesProperties.
end

