% This script performs a basic time stretcher that uses Overlap and Add
% method (OLA) without systematically correcting the phase difference 
% before the synthesis stage, but attempting to minimize the impacts of 
% phase difference for this given specific input signal, a 440Hz sine tone. 
% The reasoning for choosing suitable parameters is stated below.
% 
% Carefully choose certain N and H_A to minimize the difference between the 
% phase difference of analysis frames and the phase difference of the synthesis
% frames, i.e. we want theta_m+1[k] - theta_m[k] = phi_m+1[k] - phi_m[k] + 2*pi*m,
% m is an integer. We know that theta_m+1[k] - theta_m[k] is actually equal
% to Q*(phi_m+1[k] - phi_m[k]), hence we have phi_m+1[k] - phi_m[k] +
% 2*pi*m = Q*(phi_m+1[k] - phi_m[k]), which is (Q-1)*(phi_m+1[k] -
% phi_m[k]) = 2*pi*m. If we set N to make our f = 440Hz right on bin (e.g.,
% N = 4410, f/k = 440/44 = Fs/N = 44100/4410 = 10), then this expression will 
% be simplified to (Q-1)*H_A*w_k = 2*pi*m
% H_A = N*m/((Q-1)*k), (Q > 1 or Q < 1, N, m and k are integers)
% O = 1 - H_A/N = 1 - m/((Q-1)*k), m = 1, 2, ..., (Q-1)*k-1 
% 
% When Q = 1.25, O = 1 - m/11, if we let m = 1, then O = 10/11, H_A = N/11
% = 401
% 
% When Q = 0.75, O = 1 + m/11, if we let m = -1, then O = 10/11, H_A = N/11
% = 401
%
% This approach can also be applied to a general music file with a known 
% pitch. I tried 'Cath_cut.wav', a single melody line with the frequency of 
% the starting note as ~280H, which is also a on-bin frequency with N =
% Fs/10 = 4410 and k is 28 in this case, O = 1 - m/7 for Q = 1.25. 
% O = 1 + m/7 for Q = 0.75. In general, the larger N the
% better preservation of the audio quality from time streching. 
% O = 3/4 (H_A = N/4) is a good starting point for comparison with the
% approach mentioned above.

% Author: Zuyu Chen
% Date: Nov 8, 2024

clc
clear
close all
% 
% Set the Parameters

% Case I when Q = 1:
% w_t = 0.1; % Window Time in Seconds
% Q = 1; % Time strech factor, no expansion or compression
% O = 9/10; % Overlap Factor, such that N/H_A = 10 is an integer, 
%           % H_A = N/10 = 441, also an integer

% Read in an input WAV file and store it in the vector 
% fn = 'A440Hz.wav';
fn = 'Cath_cut.wav';
% determine the bin k from k/N = f0/Fs if set on bin
switch fn
    case 'Cath_cut.wav'
        k = 28; % The fundamental frequency of the starting note is 280 Hz
    case 'A440Hz.wav'
        k = 44; % The fundamental frequency is 440 Hz
    otherwise
        error('Unknown file name: %s', fn);
end
[x, Fs] = audioread(['audio_samples/' fn]);
% Combine stereo to mono chanel
x = sum(x,2)/2;

% Case II when Q = 0.75: 
Q = 0.75;  % Compression by a quarter of the orginal signal
w_t = 0.1; % Window Time in Seconds;
m = -1;    % chose an integer ranging from -1 to -k*(Q-1)
O = 1 - m/(k*(Q - 1)); 
% O = 3/4; % a reference for comparison

% Case III: 
% Q = 1.25;  % Expansion by a quarter of the signal
% w_t = 0.1; % Window Time in Seconds;
% m = 1;     % chose an integer ranging from 1 to k*(Q-1)
% O = 1 - m/(k*(Q - 1)); 
% % O = 3/4; % a reference for comparison

% Derive N and H_A 
N = Fs*w_t; % Frame length
tolerance = 1e-10;               % Apply small tolerance adjustments to 
                                 % handle floating-point inaccuracies 
                                 % of the cell function
H_A = ceil(N*(1-O) - tolerance); % Analysis Hopsize, use ceil to round  
                                 % to the nearest integer above
H_S = ceil(Q*H_A - tolerance);   % Synthesis Hopsize, also an integer

% Create a hann window
win = 0.5*(1 - cos(2*pi*(0:N-1)/N))'; % Periodic with 1 cycle per N points

% Zero padding the input signal for performing STFT Analysis
L = length(x); % length of x
% Zero padding at the begining 
x_padded = [zeros(N, 1); x]; 
% Solve for NF: number of frames 
NF = ceil((L+N)/H_A); 
% Solve number of zeros to padd at the end:
% L + N + num_zeros_end = H_A*(NF - 1) + N
num_zeros_end = H_A*(NF-1)-L; 
x_padded = [x_padded; zeros(num_zeros_end, 1)];

% Intilize the output vector y
y = zeros(H_S*(NF - 1) + N, 1);
% Synthesize y by adding the overlapped frames
for i = 0 : (NF - 1)  
    y(i*H_S + 1:i*H_S + N) = y(i*H_S + 1:i*H_S + N) + win.*(win.* x_padded(i*H_A + 1:i*H_A + N));
end

% For the case of Q = 1
if Q == 1

% compute the constant gain factor
    % check if N/H_A is an integer
    assert(N/H_A == floor(N/H_A));
    gain = sum(win((1:N/H_A)*H_A).^2);
    
% Plot the time domain signal: input and output
    % Plot the overlaid input and output, and the error
    subplot(2,1,1)
    plot(x_padded);
    title('Input and Output Overlaid')
    hold on
    plot(y/gain);
    subplot(2,1,2)
    plot(y/gain - x_padded)
    title('Error')
end

% Plot the spectrogram: input vs output
MA1_s2751685_Chen_myspec(x_padded, Fs, N, O) % plot the input
MA1_s2751685_Chen_myspec(y, Fs, N, O)        % plot the output
% Sound check
soundsc(y,Fs)
soundsc(x, Fs)






