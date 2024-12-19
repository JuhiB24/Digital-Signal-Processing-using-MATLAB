% EMG DATASET: Combined influence of forearm orientation and muscular contraction on EMG pattern recognition 

%% Plot the first 10 seconds and the last 10 seconds of your raw data using MATLAB’s subplot() command and label the axis.

% Initializing parameters
sampling_freq = 1000;
time = 10;
samples = sampling_freq * time;   % no. of samples

load("2014_06_20_Subject01.mat")

emg_data = data_EMG;

% Selecting the first channel and first participant
emg_channel = emg_data(:, 1, 1); % [samples x 1 x 1] 

first_10s = emg_channel(1:samples);
last_10s = emg_channel(end - samples + 1:end);

time_vector = (0:samples-1) / sampling_freq;   % x-axis

figure;

% Plotting First 10 Seconds of EMG Data
subplot(2,1,1);
plot(time_vector, first_10s);
title("First 10 Seconds of EMG Data");
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;

% Plotting Last 10 Seconds of EMG Data
subplot(2,1,2);
plot(time_vector, last_10s);
title("Last 10 Seconds of EMG Data");
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;

%% From the plots above, can you tell whether you data has some bias or offset voltage? If it doesn’t, then add a DC bias of 5 units and plot the output. Write the command to remove the DC bias and plot the result. [hint: use detrend() function]

% Check for DC bias
bias = mean(emg_channel);  
time_vector = (0:length(emg_channel)-1) / sampling_freq;

% Plot the original EMG data with the bias line
figure;
subplot(2,1,1)
plot(time_vector, emg_channel);
title('Original EMG Signal with DC Bias Line - Channel 1, Participant 1');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;
hold on;
plot(time_vector, bias * ones(size(time_vector)), '-k', 'LineWidth', 1.5); % Original DC bias line
legend('EMG Signal', 'Mean (Bias)');
hold off;

% Adding a DC bias of 5 units
emg_bias = emg_channel + 5;
new_bias = mean(emg_bias); 

% Plot the EMG signal after adding 5 units DC bias
subplot(2,1,2)
plot(time_vector, emg_bias);
title('EMG Signal with 5 Units Added DC Bias - Channel 1, Participant 1');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;
hold on;
plot(time_vector, new_bias * ones(size(time_vector)), '-k', 'LineWidth', 1.5); % New DC bias line after adding 5 units
legend('EMG Signal with 5 Units Bias', 'New Mean (Bias + 5 units)');
hold off;


% Remove DC bias using detrend
emg_det = detrend(emg_bias, 'constant');

% Plot original and detrended signal with bias/mean lines
figure;

% Original signal with bias line
subplot(2,1,1);
plot(time_vector, emg_bias);
hold on;
plot(time_vector, new_bias * ones(size(time_vector)), '--k', 'LineWidth', 1.5); 
title('Original EMG Signal with Bias Line');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
legend('Original Signal', 'Mean (Bias)');
grid on;
hold off;

% Detrended signal with zero mean line
subplot(2,1,2);
plot(time_vector, emg_det);
hold on;
plot(time_vector, zeros(size(time_vector)), '--k', 'LineWidth', 1.5); % Zero mean line
title('Detrended EMG Signal (DC Bias Removed)');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
legend('Detrended Signal', 'Mean Line ~ 0');
grid on;
hold off;

% Mean of original EMG data (bias) = 0.0044. Since, the bias is too low, I
% added a 5 unit DC bias to my EMG channel whose mean now = 5.0044. Now, using
% detrend() function on the biased channel signal the mean came out to be
% ~0. 

%% Plot a power spectrum of your signal using the periodogram() function. Label the axes.  

[PSD_period,PSD_f] = periodogram(emg_channel,[],[],1000);

figure; 
subplot(2,1,1);
stem(PSD_f,2*PSD_period);
xlabel('Analog frequency (Hz)'); 
ylabel('Power Spec. Density'); 
title('Periodogram of EMG data');
xlim([0 50]);

subplot(2,1,2);
plot(PSD_f,2*PSD_period);
xlabel('Analog frequency (Hz)'); 
ylabel('Power Spec. Density'); 
title('Periodogram of EMG data');
xlim([0 50]);

% The prominent peak at approximately 10 Hz suggests that this frequency is
% where the muscle activity is most pronounced, likely representing the 
% rhythm of hand gestures or muscle contractions. The absence of significant
% power at frequencies beyond 20 Hz suggests that the EMG signal is relatively
% low-frequency, which is consistent with typical EMG signals focused on 
% voluntary muscle contractions.

%% Now plot the power spectrum using the Welch method with different overlaps 25%, 50% 80%, 95%. Comment on the results.  

segment_length = 100; % Smaller segment length for more noticeable overlap effects

overlaps = [0.25, 0.5, 0.8, 0.95]; % Overlap percentages
figure;
for i = 1:length(overlaps)
    overlap_samples = segment_length * overlaps(i); % Calculate overlap in samples
    [PSD_welch, PSD_f] = pwelch(emg_channel, segment_length, overlap_samples, 1000, sampling_freq);
    
    subplot(2,2,i);
    plot(PSD_f, 2 * PSD_welch);
    title(['Power Spectrum with ', num2str(overlaps(i) * 100), '% Overlap']);
    xlabel('Frequency (Hz)');
    ylabel('Power/Frequency (dB/Hz)');
    grid on;
end

% Increasing the overlap smooths the power spectral estimate by increasing
% the number of overlapping segments. But in this case, despite different 
% overlap settings and lower the segment length, the overall shape of the power spectrum remains similar, 
% as the dominant frequencies do not change. As we can see maximum power is
% concentrated at the lower frequencies of the EMG signal.

%% Load the MATLAB data file “filter1.mat” and apply the impulse response variable “h” contained in that file to your signal. Plot the signal before and after filtering. What does this mystery filter do? 

% Load the filter impulse response from filter1.mat
load('filter1.mat'); % This loads the impulse response `h`

% Apply the filter to the signal using convolution
filtered_signal = conv(emg_channel, h, 'same'); 

% Time vector for plotting
time_vector = (0:length(emg_channel)-1) / sampling_freq;

% Plot the original and filtered signals
figure;

% Original Signal
subplot(2,1,1);
plot(time_vector, emg_channel);
title('Original EMG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;

% Filtered Signal
subplot(2,1,2);
plot(time_vector, filtered_signal);
title('Filtered EMG Signal');
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
grid on;

% From the plots we observe the filter seems to be a low-pass filter, 
% which attenuates high-frequency noise while preserving the 
% lower-frequency components of the EMG signal. This enhances the 
% signal quality by removing high-frequency noise or artifacts, making 
% the signal more suitable for analysis. 


%% Based on your signal type, identify the relevant signal frequency or frequency band as well as noise frequency or frequencies. What is the magnitude of the signal power and noise power? 

% Assuming the frequencies at lower frequencies below 10Hz to be the relevant frequencies. 
% Assuming the noise at around 12.5Hz (prominent peak) and above 20Hz

% From the periodogram plotted in Q6. (a), we can compute the values of signal power and noise power 
signal_band = [0, 10];       % Signal frequencies below 10 Hz
noise_band = [10, 15];       % Noise frequencies around 12.5 Hz (above 10 Hz)

% Calculate signal power by summing PSD values in the signal band
signal_indices = (PSD_f >= signal_band(1)) & (PSD_f <= signal_band(2));
signal_power = sum(PSD_period(signal_indices)) * (PSD_f(2) - PSD_f(1));

% Calculate noise power by summing PSD values in the noise band
noise_indices = (PSD_f >= noise_band(1)) & (PSD_f <= noise_band(2));
noise_power = sum(PSD_period(noise_indices)) * (PSD_f(2) - PSD_f(1));

% Display results
fprintf('Signal Power: %f\n', signal_power);
fprintf('Noise Power: %f\n', noise_power);

% Signal Power: 0.000029
% Noise Power: 0.000024

%% Design an FIR or IIR filter to remove the noise frequencies from your signal. What filter order have you selected and why? Plot the signal in time domain and its power spectrum before and after filtering.

design = filter(lpf,emg_channel);
time_vector = (0:length(emg_channel)-1) / sampling_freq;

% Plot the signal in time domain
figure;
subplot(2,1,1);
plot(time_vector, emg_channel)
xlabel('Time (seconds)');
ylabel('Amplitude (µV)');
title('EMG data before filtereing');
grid on;

subplot(2,1,2)
plot(time_vector, design)
xlabel('Time (seconds)');
ylabel('Amplitude (µV)'); 
title('EMG data after filtereing');
grid on;

% Plot power spectrum of the signal before and after filtering through Blackman LPF 

[PSD_period,PSD_f] = periodogram(emg_channel,[],[],1000);

figure; 
subplot(2,1,1);
stem(PSD_f,2*PSD_period);
xlabel('Analog frequency (Hz)'); 
ylabel('Power Spec. Density'); 
title('Periodogram of EMG data before filtereing');
xlim([0 50]);

[PSD_period,PSD_f] = periodogram(design,[],[],1000);

subplot(2,1,2);
stem(PSD_f,2*PSD_period);
xlabel('Analog frequency (Hz)'); 
ylabel('Power Spec. Density'); 
title('Periodogram of filtered EMG data');
xlim([0 50]);

% Designed an FIR, low pass filter to remove the noise frequencies from the EMG signal.
% Filter order: 1000
% A high filter order of 1000 is selected to achieve a sharp transition band between the passband and the stopband in FIR filters. Since this is a lowpass filter with a cutoff frequency of 10 Hz, the high order ensures:
% Better Frequency Resolution: The filter can more precisely distinguish between frequencies below and above the cutoff.
% Steep Roll-Off: A high-order filter achieves a steep roll-off, allowing it to attenuate frequencies above 10 Hz more effectively. This is particularly useful if you need a strict separation between the signal (below 10 Hz) and noise (above 10 Hz).
