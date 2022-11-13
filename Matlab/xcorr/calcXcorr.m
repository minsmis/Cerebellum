function xcorrValues = calcXcorr(ss1, ss2, startingTimestamp, window)

    % NOTE:
    % Calculate crosscorrelation with timestamps.
    %
    % PARAMETERS:
    % xcorrValues = Output, The results from correlation.
    % ss1 and ss2 (Array) = Arrays of spike timestamps.
    % startingTimestamps (Integer) = First timestamp of the recording.
    % window (Array) = The time window which is the correlations were calculated.

    %% Add Neuralynx MEX
    addpath('G:\ShuttleDrive\Nlx2NRD\Nlx'); % Own directory of neuralynx mex

    %% Declare variables
    samplingFrequency = 32000;

    %% Align traces
    sample1 = alignTraces(ss1, startingTimestamp, samplingFrequency);
    sample2 = alignTraces(ss2, startpoint, samplingFrequency);

    %% Calculate crosscorrelation using timestamps
    [xcorrValues, idx1, idx2] = tsXcorr(sample2, sample1, window, samplingFrequency);
end