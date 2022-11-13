function void = readLfp2Matlab(sampleChannel, refChannel)
    
    % NOTE:
    % For PSORT.
    % Export neuralynx .ncs data to .mat format with manual referencing.
    % 
    % PARAMETERS:
    % sampleChannel (Array or Number) = A channel which want to extract.
    % refChannel (Number) = A channel which want to set as reference channel (0 == pass)
    
    %% Add Neuralynx MEX
    addpath('G:\ShuttleDrive\Nlx2NRD\Nlx'); % Own directory of neuralynx mex

    %% Import LFP
    ADBitVolt = 0.000000122070312500000003;
    ADBituV = ADBitVolt * 1000000;
    fields = [1, 0, 1, 0, 1];
    extractHeader = 1;
    files=dir('*.ncs');

    channels = sampleChannel;
    references = refChannel;

    for n = 1:length(channels)

        [cscTimestamps, cscSampleFrequencies, cscSamples, cscHeader] = Nlx2MatCSC(['CSC', num2str(channels(n)), '.ncs'], fields, extractHeader, 1);
        [trialSamples(:, channels(n)), cscTimestamps] = straightenCSC(cscSamples, cscTimestamps);

        if references(n) == 0 % No references
            
            fname = strcat('csc', num2str(channels(n)), '.mat');
            ch_data = ch_data(:, channels(n)) * ADBituV;
        
        else

            fname = strcat('csc', num2str(channels(n)), '_Ref', num2str(references(n)), '.mat');
            
            [refTimestamps, refSampleFrequencies, refSamples, refHeader] = Nlx2MatCSC(['CSC', num2str(references(n)), '.ncs'], fields, extractHeader, 1);
            [refTrialSamples(:, references(n)), refTimestamps] = straightenCSC(refSamples, refTimestamps);

            ch_data = (ch_data(:, channels(n)) - refTrialSamples(:, references(n))) * ADBituV;

        end

        ch_time = [1:length(ch_data)];
        sample_rate = cscSampleFrequencies(1);
        save(fname);

    end
end