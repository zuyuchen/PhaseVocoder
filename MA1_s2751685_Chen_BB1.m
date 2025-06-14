%% 
% This script performs Phase-Vocoder-based TIME STRECHTING with a
% hann-shaped varying strech factor
% 
% Author: Zuyu Chen
% Date: Nov 24, 2024

clc
clear
close all

%% 

% Set parameters

w_t = 0.1;  % Window Time in Seconds;
O = 3/4;    % Overlap factor

% Read in an input WAV file and store it in the vector 

% fn = 'A440Hz.wav';
% fn = 'Cath_cut.wav';
% fn = 'birchcanoe.wav';
% fn = 'bodhran-cutM.wav';
% fn = 'Godin4_44.wav';
% fn = 'KS_example.wav';
fn = 'lathe.wav';
% fn = 'myks.wav';
[x, Fs] = audioread(['audio_samples/' fn]);
% Combine stereo to mono chanel
x = sum(x,2)/2;

%% 

% Other dependent parameters
N = Fs*w_t;             % Frame size
H_A = round(N*(1 - O)); % set fixed analysis hopsize

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

% Set Q as a function of time
t = linspace(0, 1, length(x_padded));  % normalized time vector from 0 to 1
Q = 0.5*(1 - 1*cos(2*pi*t))+0.5;    % hann-shaped time strecher
figure
plot(Q)
title('Hann-shaped streching factor Q(t)')
% Set the synthesis hopsize to be a vector
H_S = round(H_A*Q(H_A*(0:NF-1)+1));  % length: NF

%% 

% Initialization 
% 
NFFT = 2^nextpow2(N);
% Intilize the output vector y
y = zeros(sum(H_S) + N, 1); 

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
    theta_M = theta_m + H_S(j+1) * w_inst;

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
    y(sum(H_S(1:j+1))  + 1:sum(H_S(1:j+1)) + N) = y(sum(H_S(1:j+1)) + 1:sum(H_S(1:j+1)) + N) + win.*y_M;
    
    % resample each output frame 
    % R_A = 1/Q(j+1);
    % y_M_resampled = resample(y_M, R_A*H_S(j+1), H_S(j+1));
    % y(sum(H_S(1:j+1))  + 1:sum(H_S(1:j+1)) + N*R_A) = y(sum(H_S(1:j+1)) + 1:sum(H_S(1:j+1)) + N*R_A) + win.*y_M;
      
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
title('time strech before')
MA1_s2751685_Chen_myspec(y, Fs, N, O)        % plot the output
title('time strech after')

% Sound check
soundsc(y,Fs)
soundsc(x, Fs)
%% 

% Help Function for Phase Wrapping
% Wrap the phase to [-pi, pi] 
function wrapped_phase = ppa(phase)

    wrapped_phase = mod(phase + pi, 2*pi) - pi;

end