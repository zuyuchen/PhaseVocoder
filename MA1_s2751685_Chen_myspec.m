%% 
% This function plots the spectrogram. 
% 
% Inputs: x, a mono signal (treated as if zero padded) 
%         Fs, the sample rate; N, the frame length, 
%         O, overlap factor, a real number between 0 and 1).
% 
% Author: Zuyu Chen
% Date: Nov 24, 2024

function MA1_s2751685_Chen_myspec(x, Fs, N, O)
    
    % Check the dimension/orientation of x, convert it to a column vector
    % if it's not
    if ~iscolumn(x)
        x = x(:);
    end

    % Compute Hopsize and Number of frames
    H_A = round(N*(1-O));     % Analysis Hopsize
    NF =  ceil((length(x) - N)/H_A);        % Number of the frames
    % Create a Hann window
    win = 0.5*(1 - cos(2*pi*(0:N-1)/N))'; 

    % Initialize the first column of the STFT matrix
    NFFT = 2^nextpow2(N); 
    S = zeros(NFFT, NF);
    S(:, 1) = fft(win.* x(1:N), NFFT);
    % Complete the STFT matrix
    for i = 1 : (NF - 1)
        S(:,i+1) = fft(win.* x(i*H_A + 1:i*H_A + N), NFFT); 
    end

    %% 
    % Plot Spectrogram
    % 
    % Time axis
    tt = H_A*(0:NF-1)/Fs;  
    % Frequency axis (0 to Nyquist frequency)
    ff = linspace(0, Fs/2, NFFT/2 + 1); 
    % Convert STFT magnitude to decibels
    magnitude_dB = 20 * log10(abs(S)); 
    % discard the bins above the Nyquist bin
    magnitude_dB = magnitude_dB(1:NFFT/2+1, :); 
    
    % Plot the spectrogram
    figure;
    imagesc(tt, ff, magnitude_dB); % x-axis: time, y-axis: frequency
    axis xy;                       % Flip y-axis so it ascends upward
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');

    % Add colorbar and limit the color axis range to 60 dB updown
    caxis([max(magnitude_dB(:)) - 60, max(magnitude_dB(:))]); 
    colorbar;
    % Add the title
    text = sprintf('STFT with frames of %0.2f seconds long and %0.1f%% overlap', N/Fs, O*100);
    title(text);
end


