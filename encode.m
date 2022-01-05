function [acoeffs, gains, max_a, max_gain, pitch, max_p] = encode(...
    signal, p, fs, stepSize,...
    windowSize, numBits...
)
% encode performs LPC analysis on a signal.
% Inputs:
%   signal (vector[double]): The signal to analyze.
%   p (int): The number of LPC coeffients to generate per frame.
%   fs (int): The signal sampleing rate.
%   stepSize(int): The number of samples that seperate each frame.
%   windowSize(int): The number of samples per frame.
%   numBits(int): The number of bits used to repersent the LPC analysis values.
% Outputs:
%   acoeffs (matrix[int]): A matrix of quantized LPC coeffiencts. Each row contains the coefficents
%       for a frame
%   gains (vector[int]): A vector of quantized gains for each frame.
%   max_a (double): The maximum value of the acoeffs matrix before quantization. Used for 
%       converting the acoeffs matrix back into doubles.
%   max_gain (double): The maximum value of the gain vector before quantizaiton.
%   pitch (vector[int]): A vector of quantized pitch values.
%   max_p (double): The maximum value of the pitch vector before quantization.
    if nargin ~= 6
        numBits = 12;
    end
    % calculate teh frame overlap
    overlap = windowSize - stepSize;%
    % calcuate the pitch for each frame
    pitch = pitchdetector(signal, fs, windowSize, overlap);

    % High frequency emphasis filtering
    signal = filter([1 -0.93], 1, signal);
    % create the hamming window
    w = hamming(windowSize);
    % segment the signal into overlapping frames
    % and window each frame
    buff = w.*buffer(signal, windowSize, overlap, 'nodelay');
    % get the number of frames
    stepCount = size(buff, 2);
    % preallocate the LPC coefficent matrix and gain vector  
    acoeffs = zeros(stepCount, p);
    gains = zeros(stepCount, 1);
    for i = 1:stepCount
        % get frame out of buffer
        chunk = buff(:, i);
        % calculate lpc coffienct and residual variance
        [a, g] = lpc(chunk, p);
        % store LPC coefficents
        acoeffs(i, :) = a(2:end); % chop off first element in a, it will alway be 1
        gains(i) = sqrt(g); % gain is the std of the residual which is the square root of the variance
    end
    % extract max values for quantization
    max_a = max(abs(acoeffs(:)));
    max_gain = max(abs(gains));
    max_p = max(pitch);
    % quantize values
    acoeffs = uencode(acoeffs, numBits, max_a);
    gains = uencode(gains, numBits, max_gain);
    pitch = uencode(pitch, numBits, max_p);
end

