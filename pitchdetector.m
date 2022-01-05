function pitch = pitchdetector(...
    signal,...
    sampling_rate,...
    frame_length,...
    frame_overlap...
)
    w = hamming(frame_length);

    buff = w.*buffer(signal, frame_length, frame_overlap, 'nodelay');
    frame_count = size(buff, 2);
    pitch = zeros(frame_count, 1);

    for i = 1:frame_count
        frame = buff(:, i);
        [frame_corr, lags] = xcorr(frame);
        % isolate the cross correlation to
        % only the positive lags
        % the cross correlation is symmetric
        % around zero
        frame_corr = frame_corr(lags >= 0);
        lags = lags(lags >= 0);
        % find the peaks in the autocorrelation
        [pks, pk_loc] = findpeaks(frame_corr);
        % if there is a peak found
        % otherwise the frame is actually unvoiced
        if ~isempty(pks)
            [~, peaks_sorted_ndx] = sort(pks, 'descend');
            pitch_period = lags(pk_loc(peaks_sorted_ndx(1)));
            pitch(i) = sampling_rate/(2*pitch_period);

        end
    end
    pitch_thresh = 300;
    pitch(pitch > pitch_thresh) = 0;
end
