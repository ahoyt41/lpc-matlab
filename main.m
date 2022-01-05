clear; close all;
addpath('PESQ/obj_evaluation');

if exist('output', 'dir') == 0
    mkdir('output');
end
if exist('temp', 'dir') == 0
    mkdir('temp');
end
if exist("plots", 'dir') == 0
    mkdir ("plots");
end

% Part A)
p = 18; % number of LPC coefficients per frame
numBits = 12; % number of bits to compress LPC coeffients into
stepDuration = 0.02; % step size in seconds
windowDuration = 0.03; % window size in seconds
generateMetrics(p, numBits, stepDuration, windowDuration);


%pool = parpool('local');
%test_num_coeffs("mysample.wav");
%test_num_bits("mysample.wav");
%test_step_size("mysample.wav");
%test_window_size("mysample.wav");
%delete(pool);

p = 18;
numBits = 10;
stepDur = 0.02;
winDur = stepDur*1.25;

[mySig, fs] = audioread("mysample.wav");
stepLen = round(fs*stepDur);
winLen = round(fs*winDur);

[a, g, max_a, max_g, pitch, max_p] = encode(mySig, p, fs, stepLen, winLen, numBits);
mysig_received = decode(a, max_a, g, max_g, fs, stepLen, winLen, pitch, max_p, numBits);
mysig_received(mysig_received > 4*sqrt(var(mysig_received))) = 0;
mysig_received = mysig_received./max(abs(mysig_received));
% add aditional gain to the reconstructed signal. it is very quiet likely due to the snr
audiowrite("mysample-received.wav", mysig_received, fs);

pesq("mysample.wav", "mysample-received.wav")
total_bits = numBits*(numel(a) + numel(g) + numel(pitch)) + 192;
bits_per_sample = total_bits / length(mySig) 


function test_num_bits(inputfile)
    p = 10;
    [x, fs] = audioread(inputfile);
    stepSize = round(0.02*fs);
    windowSize = round(0.03*fs);
    num_bits_vector = 8:16;
    pesq_vector = zeros(size(num_bits_vector));
    bits_per_sample_vector = zeros(size(num_bits_vector));
    parfor i = 1:length(num_bits_vector)
        num_bits = num_bits_vector(i);
        [a, g, max_a, max_g, pitch, max_p] = encode(...
            x, p, fs, stepSize, windowSize, num_bits);
        y = decode(...
            a, max_a, g, max_g, fs, stepSize, ...
            windowSize, pitch, max_p, num_bits);
        outputfile = sprintf("temp/recreated-%d.wav", i);
        audiowrite(outputfile, y, fs);
        pesq_vector(i) = pesq(inputfile, outputfile);
        total_bits = num_bits*(numel(a) + numel(g) + numel(pitch)) + ...
            192;
        bits_per_sample_vector(i) = total_bits/length(x);
    end
    figure
    plot(num_bits_vector, pesq_vector);
    xlabel("Bits per Coefficent");
    ylabel("PESQ");
    hold on
    yyaxis right;
    plot(num_bits_vector, bits_per_sample_vector);
    ylabel("Bits per Sample");
    hold off
    title(sprintf("PESQ and Bits per Sample with\nVarying Quantization Levels"));
    saveas(gcf, "plots/quant-levels.jpg");
    delete("temp/*.wav");
end

function test_num_coeffs(inputfile)
    p_vector = 8:64;
    pesq_vector = zeros(size(p_vector));
    bits_per_sample_vector = zeros(size(p_vector));
    [x, fs] = audioread(inputfile);
    stepSize = round(0.02*fs);
    windowSize = round(0.03*fs);
    num_bits = 12;
    parfor i = 1:length(p_vector)
        p = p_vector(i);
        [a, g, max_a, max_g, pitch, max_p] = encode(...
            x, p, fs, stepSize, windowSize, num_bits);
        y = decode(...
            a, max_a, g, max_g, fs, stepSize, ...
            windowSize, pitch, max_p, num_bits);
        outputfile = sprintf("temp/recreated-%d.wav", i);
        audiowrite(outputfile, y, fs);
        pesq_vector(i) = pesq(inputfile, outputfile);
        total_bits = num_bits*(numel(a) + numel(g) + numel(pitch)) + ...
            192;
        bits_per_sample_vector(i) = total_bits/length(x);
    end
    figure
    plot(p_vector, pesq_vector);
    xlabel("Number of LPC Coeffients");
    ylabel("PESQ");
    hold on
    yyaxis right;
    plot(p_vector, bits_per_sample_vector);
    ylabel("Bits per Sample");
    hold off
    title(sprintf("PESQ and Bits per Second\nFor Diffferent p Values"));
    saveas(gcf, "plots/p-vals.jpg")
    delete("temp/*.wav");
end

function test_step_size(inputfile)
    p = 10;
    [x, fs] = audioread(inputfile);
    stepSizeVector = round((0.01:0.005:0.03)*fs);
    windowSizeVector = round(1.25*stepSizeVector);
    num_bits = 12;
    pesq_vector = zeros(size(stepSizeVector));
    bits_per_sample_vector = zeros(size(stepSizeVector));

    parfor i = 1:length(stepSizeVector)
        stepSize = stepSizeVector(i);
        windowSize = windowSizeVector(i);
        [a, g, max_a, max_g, pitch, max_p] = encode(...
            x, p, fs, stepSize, windowSize, num_bits);
        y = decode(...
            a, max_a, g, max_g, fs, stepSize, ...
            windowSize, pitch, max_p, num_bits);
        outputfile = sprintf("temp/recreated-%d.wav", i);
        audiowrite(outputfile, y, fs);
        pesq_vector(i) = pesq(inputfile, outputfile);
        delete(outputfile);
        total_bits = num_bits*(numel(a) + numel(g) + numel(pitch)) + ...
            192;
        bits_per_sample_vector(i) = total_bits/length(x);
    end
    figure
    plot(stepSizeVector/fs, pesq_vector);
    xlabel("Step Size (s)");
    ylabel("PESQ");
    hold on
    yyaxis right;
    plot(stepSizeVector/fs, bits_per_sample_vector);
    ylabel("Bits per Sample");
    hold off
    title(sprintf("PESQ and Bits per Sample\nfor different Step Sizes"));
    saveas(gcf, "plots/step-size.jpg");
    delete("temp/*.wav");
end

function test_window_size(inputfile)
    [x, fs] = audioread(inputfile);
    p = 10;
    stepSize = round(0.02*fs);
    windowSizeVector = round((1:0.01:2)*stepSize);
    num_bits = 12;
    pesq_vector = zeros(size(windowSizeVector));
    bits_per_sample_vector = zeros(size(windowSizeVector));

    parfor i = 1:length(windowSizeVector)
        windowSize = windowSizeVector(i);
        [a, g, max_a, max_g, pitch, max_p] = encode(...
            x, p, fs, stepSize, windowSize, num_bits);
        y = decode(...
            a, max_a, g, max_g, fs, stepSize, ...
            windowSize, pitch, max_p, num_bits);
        outputfile = sprintf("temp/recreated-%d.wav", i);
        audiowrite(outputfile, y, fs);
        pesq_vector(i) = pesq(inputfile, outputfile);
        total_bits = num_bits*(numel(a) + numel(g) + numel(pitch)) + ...
            192;
        bits_per_sample_vector(i) = total_bits/length(x);
    end
    figure
    plot(windowSizeVector/fs, pesq_vector);
    xlabel("Window Duration (s)");
    ylabel("PESQ");
    hold on
    yyaxis right;
    plot(windowSizeVector/fs, bits_per_sample_vector);
    hold off
    title(sprintf("PESQ and Bits per Sample for differing Frame Sizes"));
    saveas(gcf, "plots/frame-size.jpg");
    delete("temp/*.wav");
end

function metrics = generateMetrics(p, numBits, stepDuration, windowDuration)
    % generatMetrics analyzes and resynthesizes all of the wave files found in the Project2/
    % directory. It returns a table with metrics for each wave file and saves
    % it to a file name metrics.csv.
    % Inputs:
    %   p (int): The number of LPC coeffients to generate per frame
    %   numBits (int): The number of bits each LPC coeffiecient is represented with
    %   stepDuration (double): The number of second that each frame is seperated by
    %   windowDuration (double): The length of an analysis frame in seconds.
    % Output:
    %   metrics (table): A table containing the metrics for a a file.
    %   talbe columns:
    %       Reconstructe File Name
    %       Bits: The total number of bits used for the LPC analysis
    %       Samples per second: The number of samples needed to create one second of audio
    %       Bits per second: The number of bits needed to create one second of audio
    %       Compression ration: The ratio between the number of bits in the LPC analysis
    %           and the number of bits in the original file.
    %       PESQ: The PESQ score of the reconstructed signal.

    % fetch all of the audio files from the project 2 directory
    inputFiles = dir("Project2/*.wav");
    % initalize metric vectors
    ReconstructedFile = cell(length(inputFiles), 1);
    Bits = zeros(length(inputFiles), 1);
    SamplesPerSecond = zeros(length(inputFiles),1);
    BitsPerSecond = zeros(length(inputFiles), 1);
    PESQ = zeros(length(inputFiles), 1);
    CompressionRatio = zeros(length(inputFiles), 1);
    for i = 1:length(inputFiles) % for each file
        % get the inputfile name and create the outputfile name
        % as output/Sample<number>-receieved.wav
        filename = inputFiles(i).name;
        name_no_ext = strsplit(filename, '.');
        inputfile = fullfile("Project2", filename);
        outputfile = fullfile("output", sprintf("%s-received.wav", name_no_ext{1}));

        % extract the file metadata for the input fie
        % and read the data for the input file
        inputInfo = audioinfo(inputfile);
        [inputSignal, fs] = audioread(inputfile);

        % calcuate the number of sample per step and per frame
        % from the duration and sample rate
        stepSize = round(stepDuration*fs);
        windowSize = round(windowDuration*fs);

        % calculate the LPC coeffiencts, gain, and pitch values
        [a, gain, max_a, max_g, pitch, max_p] = encode(inputSignal, p, fs, stepSize, windowSize, numBits);
        % synthesize the reconstructed signal from the LPC coeffiencts
        outputSignal = decode(a, max_a, gain, max_g, fs, stepSize, windowSize, pitch, max_p, numBits);
        % save the reconstructed file
        audiowrite(outputfile, outputSignal, fs);
        % Calculate the metrics for the LPC analysis and reconstruction
        ReconstructedFile{i} = outputfile;
        % Each LPC coeffient, gain, and pitch values is represented
        % by numBits. We also need three double values to store the maximum
        % coffient, gain, and pitch value for ideal reconstruction which. These
        % three doubles take up 192 bit, 64 each
        Bits(i) = numBits*(numel(a) + numel(gain) + numel(pitch)) + 192;
        SamplesPerSecond(i) = (numel(a) + numel(gain) + numel(pitch)) * (length(inputSignal)/fs);
        bits_per_sample = Bits(i)/length(inputSignal);
        BitsPerSecond(i) = Bits(i)/(length(inputSignal)/fs);
        PESQ(i) = pesq(inputfile, outputfile);
        CompressionRatio(i) = bits_per_sample/inputInfo.BitsPerSample;
    end
    % create the table, display it, and write it to a file.
    metrics = table(ReconstructedFile, Bits, SamplesPerSecond, BitsPerSecond, CompressionRatio, PESQ);
    disp(metrics);
    writetable(metrics, 'metrics.csv');
end
