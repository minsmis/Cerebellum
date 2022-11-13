function [tsOffsets, ts1idx, ts2idx] = crosscorrelogram(ts1, ts2, window)
    % INPUT:
    % ts1 - timestamp array, center of cross-correlogram
    % ts2 - timestamp array
    % window - the window of the timestamps to compute, e.g. [-0.1 0.1] for
    %           100 milliseconds around each spike.
    %
    % OUTPUT:
    % tsOffsets - the offsets from each spike in ts1 that has spikes nearby
    % ts1idx - the index into ts1 for each offset in tsOffsets
    % ts2idx - the index into ts2 for each offset in tsOffsets
	%
	% If you want to make a cross-correlogram histogram, you can make a bar plot
	% hist(tsOffsets, 100);
	%
	% NOTE:
	% If you don't need millisecond precision, you're probably better off 
	% binning your spikes and using xcorr. If you're hardcore (you're hardcore
	% aren't you?) then use this function.
    %
    % LICENSE:
    % Copyright (c) 2010, Alex Wiltschko
    % All rights reserved.
    %
    % Redistribution and use in source and binary forms, with or without
    % modification, are permitted provided that the following conditions are
    % met:
    %
    %   * Redistributions of source code must retain the above copyright
    %     notice, this list of conditions and the following disclaimer.
    %   * Redistributions in binary form must reproduce the above copyright
    %     notice, this list of conditions and the following disclaimer in
    %     the documentation and/or other materials provided with the distribution
    %
    % THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
    % AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
    % IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
    % ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
    % LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
    % CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
    % SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
    % INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
    % CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
    % ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
    % POSSIBILITY OF SUCH DAMAGE.
	
	% Make sure that all arrays are single rows, not columns
	if size(ts1,1) > size(ts1,2); ts1 = ts1'; end
	if size(ts2,1) > size(ts2,2); ts2 = ts2'; end
		
    startWindow = 1;
    spikeIdx = 2;
    ccInfo = zeros(length(ts1), 3)+NaN; % startIdx, endIdx, spikeTime, instantRate
	dumb = 0;
	
    while spikeIdx <= length(ts1)

        % Seek to the beginning of the current spike's window
        i = startWindow;
        while ts2(i) <= (ts1(spikeIdx) + window(1)) && i < length(ts2)
            i = i+1;
        end;
		
		startWindow = i; % save the location for later

        if ts2(i) > ( ts1(spikeIdx) + window(2))
            spikeIdx = spikeIdx+1;
            continue;
        end
        
        % Find all the spike indices that fall within the window
        while ts2(i) <= (ts1(spikeIdx) + window(2)) && i < length(ts2)
            i = i+1;
        end
        endWindow = i-1;
        
        ccInfo(spikeIdx,1) = startWindow;
        ccInfo(spikeIdx,2) = endWindow;
        ccInfo(spikeIdx,3) = ts1(spikeIdx);
        spikeIdx = spikeIdx+1;
		
    end
        
    % Now I've got all of the indices into spikes, and their offset
    % from the center-spike is easily calculable.
    % I think I'll increment into the array bit by bit.
    diffArray = diff(ccInfo (:,1:2),1,2); 
    done = 0; tsOffsets = []; instantRate = []; ts1idx = []; ts2idx = [];
    incr = 0;
    
    while ~done
        tmp = find(diffArray >= incr);
        if isempty(tmp); done = 1; continue; end
        idx = ccInfo(tmp,1)+incr; % brilliant!
        centerTimes = ccInfo(tmp,3);
        tsOffsets = [tsOffsets ts2(idx)-centerTimes'];
		
		ts1idx = [ts1idx; tmp];
		ts2idx = [ts2idx; idx];
        incr = incr+1;
    end
 
end