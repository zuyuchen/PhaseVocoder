%% 
% This script performs time streching with Phase Vocoder 
% 
% Author: Zuyu Chen
% Date: Nov 14, 2024

clc
clear
close all

%% 

% Case I when Q = 1:
% 
% w_t = 0.1; % Window Time in Seconds
% Q = 1; % Time strech factor, no expansion or compression
% O = 9/10; % Overlap Factor, such that N/H_A = 10 is an integer, 
% %           % H_A = N/10 = 441, also an integer
%% 

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

%% 

% Case II when Q = 0.75: 

Q = 0.75;  % Compression by a quarter of the orginal signal
w_t = 0.1; % Window Time in Seconds;
m = -1;    % chose an integer ranging from -1 to -k*(Q-1)
O = 1 - m/(k*(Q - 1)); 
% O = 3/4; % a reference for comparison
%% 
% 
% Case III when Q = 1.25: 

% Q = 1.25;  % Expansion by a quarter of the signal
% w_t = 0.1; % Window Time in Seconds;
% m = 1;     % chose an integer ranging from 1 to k*(Q-1)
% O = 1 - m/(k*(Q - 1)); 
% % O = 3/4; % a reference for comparison

%% 
% 
% Case IV when Q = 5, 
% H_S = 5*H_A = 5*(N*(1-O)) < N, O needs to be larger than 4/5, 
% H_A needs to be no greater than N/5

% Q = 5;  % Expansion by a quarter of the signal
% w_t = 0.1; % Window Time in Seconds;
% m = 1;     % chose an integer ranging from 1 to k*(Q-1)
% O = 1 - m/(k*(Q - 1)); % O needs to be at least larger than 4/5 
%% 
% 
% Derive N H_A 
N = Fs*w_t; % Frame length
tolerance = 1e-10;               % Apply small tolerance adjustments to 
                                 % handle floating-point inaccuracies 
                                 % of the cell function
H_A = ceil(N*(1-O) - tolerance); % Analysis Hopsize, use ceil to round  
                                 % to the nearest integer above
H_S = ceil(Q*H_A - tolerance);   % Synthesis Hopsize, also an integer
%% 
% 
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

%% 

% Initialization 

NFFT = 2^nextpow2(N);
% Intilize the output vector y
y = zeros(H_S*(NF - 1) + N, 1);
% Initialize the first frame
j = 0;
x_m = win.*(x_padded(j*H_A + 1:j*H_A + N));
X_m = fft(x_m, NFFT); % should be all zeros because of the padding
% Extract the phase of the first frame
phi_m = angle(X_m);          % also should be zeros 
phi_m = phi_m(1:NFFT/2+1);   % keep only the first 1 + NFFT/2 bins
theta_m = phi_m;             % modified phase (no modification on the first frame)
% Define the nominal angular frequency w_k (only need NFFT/2+1 DFT bins)
k = (0:NFFT/2)'; 
w_k = 2*pi*k/NFFT; % 0 <= w_k <= pi

% Loop over the remaining frames
for j = 1 : (NF - 1)  

    % Perform DFT on each analysis frame
    X_M = fft(win.*x_padded(j*H_A + 1:j*H_A + N), NFFT);
%% 

    % disp('DC component of X_M')
    % disp(X_M(1));
    % disp('Nyquist frequency bin of X_M')
    % disp(X_M(4097));
%% 

    % Separate magnitude and phase
    Xmag = abs(X_M);
    Xmag = Xmag(1:NFFT/2 + 1); % keep only the first NFFT/2 + 1 bins
    Xang = angle(X_M);
    phi_M = Xang(1:NFFT/2 + 1); % keep only the first NFFT/2 + 1 bins

    % Estimate the instantaneous frequency
    w_inst = w_k + ppa(phi_M - phi_m - H_A*w_k)/H_A;
    % Calculate the modified phases for the current frame
    theta_M = theta_m + H_S * w_inst;

    % Create the new modified DFT frame of full length NFFT
    % The first half bins k = 0 ... NFFT/2
    Y_M_left = Xmag(1:end-1).* exp(1i*theta_M(1:end-1));  
    % The remaining bins k = NFFT/2+1 ... NFFT-1 are hermitian symmetric
    Y_M_right = Xmag(end-1:-1:2).* exp(-1i*theta_M(end-1:-1:2));
    % Concatenate together with the nyquist bin
    Y_M = [Y_M_left; real(Xmag(end)*exp(1i*theta_M(end))); Y_M_right];
    
    % disp('DC component of Y_M')
    % disp(Y_M(1));

    % Ensure the DC component after modification remains real
    if ~isreal(Y_M(1))
        Y_M(1) = real(Y_M(1));
    end
    % disp('Nyquist frequency component of Y_M')
    % disp(Y_M(NFFT/2 + 1));
  

    % IDFT
    y_M = ifft(Y_M, NFFT);
    assert(isreal(y_M), 'The frame of ifft is not real')
    % Synthesize y by adding the overlapped frames
    y_M = y_M(1:N); % Truncate the frame to N point
    y(j*H_S + 1:j*H_S + N) = y(j*H_S + 1:j*H_S + N) + win.*y_M;
    
    % update the phases
    theta_m = theta_M;
    phi_m = phi_M;
    clc
end
%% 
% 
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
%% 

% Help Function for Phase Wrapping
% Wrap the phase to [-pi, pi] 
function wrapped_phase = ppa(phase)

    wrapped_phase = mod(phase + pi, 2*pi) - pi;

end