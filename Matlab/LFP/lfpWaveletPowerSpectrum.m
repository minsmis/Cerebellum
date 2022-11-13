function void = lfpWaveletPowerSpectrum(givenChannels)
    
    % NOTE:
    % This function returns LFP power and compatible frequencies which were calculated by wavelet methods
    % as .mat format (matlab varialbe format).
    %
    % PARAMETERS:
    % givenChannels = array of channels which want to get power spectrum.
    
    %% Add Neuralynx MEX
    addpath('G:\ShuttleDrive\Nlx2NRD\Nlx'); % Own directory of neuralynx mex

    %% Import LFP
    fields = [1,0,1,0,1];
    extractHeader = 1;
    files = dir('*.ncs');
    channels = givenChannels; % For 8 tetrodes
    % channels = 1:2; % Number of brain structures were recorded

    for n = 1:length(channels)

        [cscTimestamps, cscSampleFrequencies, cscSamples, cscHeader] = Nlx2MatCSC(['CSC',num2str(channels(n)),'.ncs'], fields, extractHeader, 1);
        [trialSamples(:,channels(n)), cscTimestamps] = straightenCSC(cscSamples, cscTimestamps);

    end

    %% Downsampling
    trialSamples1 = downsample(trialSamples, 16); % fs = 32000 hz --> fs = 2000 hz

    %% Power calculation using wavelet
    for n = 1:length(channels)
    
        [lfpWave, freq] = cwt(trialSamples1(:, channels(n)), 2000); % 2000 = Downsampled sampling frequency
        lfpPower(:,n) = mean((abs(lfpWave).^2),2); % 2nd dimention is power, 1st dimention is time;

    end

    %% Save outputs
    save('lfpPower.mat');
end