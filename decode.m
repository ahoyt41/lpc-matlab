function signal = decode(...
    acoeffs,...
    max_a,...
    gains,...
    max_gain,...
    fs,...
    stepSize,...
    windowSize,...
    pitch,...
    max_p,...
    numBits...
)
% decode synthesizes a signal from the LPC analysis values.
% Inputs:
%   acoeffs (matrix[int]): A matrix of qunatized LPC coeffients. Each row is the values for a frame.
%   max_a (double): The maximum value of the LPC coeffiencts before quantization.
%   gains (vector[int]): A vector of quantized gain values.
%   max_gain (double): The maximum gain value before qunatization.
%   fs (int): The signal sampling rate.
%   stepSize (int): The number of samples between frames.
%   windowSize (int): The number of samples in a frame.
%   pitch (int): A vector of quantized pitch values
%   max_p (double): The maximum pitch value before quantization
%   numBits (int): The number of bits used to quantize the LPC analysis values
% Output:
%   signal (vector[double]): The sysnthesized audio sequence.
    if nargin == 9
        numBits = 12;
    end
    
    % calculate the frame overlap
    overlap = windowSize - stepSize;
    stepCount = size(acoeffs, 1);
    % generate a hamming window for overlap add
    window = hamming(windowSize);
    % initialize the reconstructed signal vecotr
    signal = zeros(stepCount*stepSize + overlap, 1);

    % decompress LPC coeffiecents
    acoeffs = udecode(acoeffs, numBits, max_a);
    gains = udecode(gains, numBits, max_gain);
    pitch = udecode(pitch, numBits, max_p);

    for i = 1:stepCount
        % decode the given frame given the LPC coeffients, gain, and pitch
        chunk = decodeFrame(acoeffs(i,:), gains(i), fs, pitch(i), window);
        % calculate the index range for the frame
        sig_ndx = (i-1)*stepSize + 1:i*stepSize + overlap;
        % overlap add
        signal(sig_ndx) = signal(sig_ndx) + chunk;
    end
    % normalize the reconstructed signal to avoid clipping
    signal = signal./max(abs(signal));
end

function chunk = decodeFrame(acoeff, gain, fs, pitch, window)
   % decodeFrame synthesizes an audio frame given its LPC coeffients
   % Inputs:
   %    acoeff (vector[double]): The LPC coeffiencts for the frame
   %    gain (double): The gain of the frame
   %    fs (int): The sequence sampling rate
   %    pitch (double): The fundemental frequency for the frame, is 0 for unvoiced
   %    window (vector[double]): The analysis window
   % Output:
   %    chunk (vector[double]): The synthesized frame.
    
   % generate the excitation signal
    input = rand(size(window));
    % if the frame is voiced
    if pitch ~= 0
        % calulate the sample period of F0
        period = round(fs/(2*pitch));
        % add impulses at the sample period
        input(1:period:end) = input(1:period:end) + 1;
    end
    % generate the frame by filtering the excitation signal
    % with the gain and LPC coefficents
    % and window the output.
    chunk = window.*filter(gain, [1, acoeff], input);
end

