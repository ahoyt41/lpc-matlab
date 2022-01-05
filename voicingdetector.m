function voicings = voicingdetector(...
    signal,...
    sampling_rate,...
    frame_length,...
    frame_overlap...
)
    w = hamming(frame_length);
    buff = w.*buffer(signal,...
        frame_length,...
        frame_overlap,...
        'nodelay'...
    );
    num_frames = size(buff, 2);
    zcr_rate = zeros(num_frames, 1);
    energy = zeros(num_frames, 1);
    voicings = zeros(num_frames,1);

    for i = 1:num_frames
        frame = buff(:, i);
        energy(i) = var(frame);
        zcr = sum(abs(diff(sign(frame))./2));
        zcr_rate(i) = zcr*sampling_rate/frame_length;
    end

    % normalize the energy
    energy = energy./max(abs(energy));

    % set the zero crossing threshold at 10%
    % any zero crossing rating below this value
    % with sufficient energy is classified
    % as voiced
    zcr_thresh = 0.1*max(zcr_rate);
    % set the energy threshhold at 1%
    % any frame with energy below this threshold
    % is classified as silence
    energy_thresh = 0.01;

    % classify the voiced and unvoiced frames via
    % logical indexing
    voicings((zcr_rate < zcr_thresh) & (energy > energy_thresh)) = 1;
    voicings((zcr_rate > zcr_thresh) & (energy > energy_thresh)) = -1;
end
