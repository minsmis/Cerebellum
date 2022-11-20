function Q = standardxcov (sample1, sample2, histogram, bin, duration)
    
    % INTRODUCTION:
    % Calculate standardized cross-covariance from cross-correlogram.
    %
    % VARIABLES:
    % sample1: timestamps from first unit.
    % sample2: timestamps from second unit.
    % histogram: histogram of cross-correlation of two units were given.
    % bin (ms): Duration of bin used for calculating the histogram was given.
    % duration (ms): Total duration of recordings.
    %
    % NOTES:
    % The level of significance (p<0.05) was thus divided by 60 yielding a
    % critical z-value of 3.34.
    %
    % REFERENCES:
    % 1. C. de Solages et al., 2008, Neuron.
    % 2. Siapas et al., 2005, Neuron.

    n1 = length(sample1);
    n2 = length(sample2);
    
    A = (n1 * n2 * bin)/duration;

    Q = (histogram - A)/sqrt(A);

end