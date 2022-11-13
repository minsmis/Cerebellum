function [xcorrValues_ms, idx1, idx2] = tsXcorr (ts1, ts2, window, fs)
    
    % CROSS-CORRELATION WITH TIMESTAMPS
    %
    % VARIABLES: 
    % ts1 = timestamp 1-D array
    % ts2 = timestamp 1-D array
    % window = 1-D array
    % fs = integer
    %
    % USAGE GUIDE:
    % Drag and drop timestamp file (.csv or .mat).
    % Run the code.
    % Import 'crosscorr.mat' to Graphpad Prism.
    % Analyze 'frequency distribution' with 1/10 (or 1/100) bin number.
    
    % CALCULATE CROSS-CORRELOGRAM
    [xcorrValues, idx1, idx2] = crosscorrelogram(ts1, ts2, window);

    % SHOW HISTOGRAM
    figure();
    histogram(xcorrValues, fs/100); % fs/100 is arbitrary value.

    % EXPORT CROSS-CORRELOGRAM
    xcorrValues_ms = xcorrValues(:) * 1000; % Convert second scale to millisecond scale
    save('crosscorr.mat');
end