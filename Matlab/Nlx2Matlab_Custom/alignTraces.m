function samples = samplingTrace (TimeStamps,startPoint,fs)
    
    % INTRODUCTION:
    % Convert PSORT timestamp files to second scale
    % 
    % FORMULA:
    % {(each timestamp) - (timestamp of initiation)} / sampling frequency
    %
    % VARIABLES:
    % samples = Output, timestamp by second scale
    % TimeStamps = Neuralynx or Psort timestamps, int64 format
    % startPoint = (OPTION) timestamp of recording intiation
    % fs = sampling frequency
    %
    % EXAMPLES:
    % If startPoint is available --> samples = samplingTrace (TimeStamps, startPoint, fs);
    % If startPoint is not available --> samples = samplingTrace (TimeStamps, [], fs);

    if startPoint ~= []
        samples = (TimeStamps(:)-startPoint)/fs;
    else
        samples = (TimeStamps(:)-TimeStamps(1))/fs;
    end
end