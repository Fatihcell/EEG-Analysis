%% ERP Analysis of O1 Channel  Based on Visual Cue


% Generate EEG data
fs = 500; % Sampling frequency
t = EEG.times; % Time vector

% Generate synthetic EEG data by adding noise
eegData = squeeze(EEG.data(13,:,:));

% Compute the ERP by averaging across trials
erp = mean(eegData, 2);

% Plot the synthetic EEG data and ERP
figure;
subplot(2,1,1);
plot(t, eegData);
title('EEG Data');
xlabel('Time (s)');
ylabel('Amplitude');

subplot(2,1,2);
plot(t, erp, 'LineWidth', 2);
title('Averaged ERP of O1 Channel');
xlabel('Time (s)');
ylabel('Amplitude');

%% Sub-band ERPs

fs = 500; % Sampling frequency
t = EEG.times; % Time vector
nTrials = size(EEG.data, 3); % Number of trials

% Select inverval for ERP sub-band Analysis
winStart = -500;
winEnd = 500;
timeIdx = find(t >= winStart & t <= winEnd);
t       = t(timeIdx);
% Generate synthetic EEG data by extracting a channel (e.g., O1)
eegData = double(squeeze(EEG.data(13,:,:)));

% Compute the ERP by averaging across trials
erp = mean(eegData, 2);

% Define frequency bands
bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
bandLimits = [0.5 4; 4 8; 8 13; 13 30; 30 100];

% Initialize sub-band ERPs
subBandERPs = zeros(size(eegData, 1), length(bands));

% Filter data and compute sub-band ERPs
for i = 1:length(bands)
    band = bandLimits(i,:);
    % Bandpass filter the data
    [b, a] = butter(3, band/(fs/2));
    filteredData = filtfilt(b, a, eegData);
    
    % Compute the ERP for the current band
    subBandERPs(:, i) = mean(filteredData, 2);
end

% Plot the original EEG data and ERP
figure;
subplot(4,2,1);
plot(t, eegData(timeIdx,:));
title('EEG Data');
xlabel('Time (ms)');
ylabel('Amplitude');

subplot(4,2,2);
plot(t, erp(timeIdx), 'LineWidth', 2);
title('Averaged ERP of O1 Channel');
xlabel('Time (ms)');
ylabel('Amplitude');

% Plot the sub-band ERPs
for i = 1:length(bands)
    subplot(4,2,i+2);
    plot(t, subBandERPs(timeIdx,i), 'LineWidth', 2);
    title(['ERP in ', bands{i}, ' Band']);
    xlabel('Time (ms)');
    ylabel('Amplitude');
end

% Adjust the layout for better visualization
title('EEG Data and Sub-Band ERPs');
%% Cross - Correlation Analysis
%%
% Define time window for cross-correlation (-500 ms to 500 ms)
winStart = -500;
winEnd = 500;
timeIdx = find(t >= winStart & t <= winEnd);

% Initialize table to store results
results = table('Size', [nchoosek(length(bands), 2), 3], ...
    'VariableTypes', {'string', 'double', 'double'}, ...
    'VariableNames', {'Comparison', 'MaxCrossCorr', 'Lag'});

% Compare sub-band ERPs
row = 1;
figure;
for i = 1:length(bands)-1
    for j = i+1:length(bands)
        % Extract the relevant time window
        erp1 = subBandERPs(timeIdx, i);
        erp2 = subBandERPs(timeIdx, j);
        
        % Compute cross-correlation
        [crossCorr, lags] = xcorr(erp1, erp2, 'coeff');
        
        % Find maximum cross-correlation and corresponding lag
        [maxCorr, maxIdx] = max(crossCorr);
        maxLag = lags(maxIdx) * (1000 / fs); % Convert to ms
        
        % Store results in table
        results.Comparison(row) = strcat(bands{i}, ' vs ', bands{j});
        results.MaxCrossCorr(row) = maxCorr;
        results.Lag(row) = maxLag;
        
        % Plot cross-correlation
        subplot(5, 2, row);
        plot(lags * (1000 / fs), crossCorr, 'LineWidth', 2);
        title(['Cross-correlation: ', bands{i}, ' vs ', bands{j}]);
        xlabel('Lag (ms)');
        ylabel('Cross-correlation');
        
        row = row + 1;
    end
end

% Display the results table
disp(results);
%% Visualize Max. Cross- Correlation and Corresponding Time-Lags
%%

% Extract bands and initialize matrices for correlation and lag
bands = {'Delta', 'Theta', 'Alpha', 'Beta', 'Gamma'};
nBands = length(bands);
corrMatrix = zeros(nBands);
lagMatrix = zeros(nBands);

% Populate matrices with max cross-correlation values and lag values
for i = 1:nBands-1
    for j = i+1:nBands
        % Construct comparison strings
        comparison1 = [bands{i} ' vs' bands{j}];
        comparison2 = [bands{j} ' vs' bands{i}];
        
        % Find indices in results for both possible comparisons
        idx1 = find(strcmp(results.Comparison, comparison1));
        idx2 = find(strcmp(results.Comparison, comparison2));
        
        % Use the first non-empty index found
        if ~isempty(idx1)
            idx = idx1;
        elseif ~isempty(idx2)
            idx = idx2;
        else
            idx = []; % No match found
        end
        
        % Populate matrices if index found
        if ~isempty(idx)
            corrMatrix(i, j) = results.MaxCrossCorr(idx);
            corrMatrix(j, i) = results.MaxCrossCorr(idx); % Symmetric matrix
            lagMatrix(i, j) = results.Lag(idx);
            lagMatrix(j, i) = results.Lag(idx); % Symmetric matrix for lag
        end
    end
end

% Create heatmap for cross-correlation values
figure;
subplot(1, 2, 1);
heatmap(bands, bands, corrMatrix, 'Colormap', jet, 'ColorScaling', 'scaled', ...
    'FontSize', 10, 'ColorLimits', [min(results.MaxCrossCorr) max(results.MaxCrossCorr)]);
title('Max Cross-Correlation between EEG Frequency Bands');
xlabel('Frequency Bands');
ylabel('Frequency Bands');
colorbar;

% Create heatmap for lag values
subplot(1, 2, 2);
heatmap(bands, bands, lagMatrix, 'Colormap', parula, 'ColorScaling', 'scaled', ...
    'FontSize', 10, 'ColorLimits', [min(results.Lag) max(results.Lag)]);
title('Lag between EEG Frequency Bands');
xlabel('Frequency Bands');
ylabel('Frequency Bands');
colorbar;

title('Cross-Correlation and Lag between EEG Frequency Bands');