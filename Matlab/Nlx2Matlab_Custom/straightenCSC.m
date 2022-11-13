function [cscOut,ts] = straightenCSC(cscIn, cscTs)
%STRAIGHTENCSC assignes ts values to each csc data point
%  usage: [cscOut,ts] = straightenCSC(cscIn, cscTs)
% 
% input a CSC matrix (from neuralynx) and its timeStamp vector
% this program will assign a timestamp to each individual sample
%

% search for discontinuity
ind = find((diff(cscTs))>1.2*mean(diff(cscTs)));

numSamples = size(cscIn,1);
ts = zeros(size(cscIn));

if ~ind % if there are no discontinuites
    for i = 1:size(cscIn,2)
        if i<size(cscIn,2) % if this is not the last i; otherwise, use previous tsPerSample
            timeStampsPerSample = (cscTs(i+1)-cscTs(i))/numSamples;
        end
        for j = 1:numSamples % now give each sample a timestamp
            ts(j,i) = cscTs(i) + timeStampsPerSample*(j-1);
        end
    end
else  % if there are discontinuities
    for i = 1:size(cscIn,2)
        if i<size(cscIn,2) % if this is not the last i, calculate ts per Sample
                           % otherwise, use previous tsPerSample
            if ~any(ind == i) % check if this i is at a discontinuity
                              % if not, use current block to calculate ts per sample
                timeStampsPerSample = (cscTs(i+1)-cscTs(i))/numSamples;
            else %if this i is discontinuity
                k = i+1;
                while find(ind == k) % find the next that is not a discontinuity
                    k = k+1;
                end 
                % and use it to calculate ts per sample
                timeStampsPerSample = (cscTs(k+1)-cscTs(k))/numSamples;
            end
        end
        for j = 1:numSamples % now give each sample a timestamp
            ts(j,i) = cscTs(i) + timeStampsPerSample*(j-1);
        end
    end
end

cscOut = cscIn(:); % straighten the csc
ts = ts(:);        % and the corresponding timestamps

if size(cscOut,1) > 1
    cscOut = cscOut';
end
if size(ts,1) > 1
    ts = ts';
end

