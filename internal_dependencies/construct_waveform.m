% freq_i [Hz]
% fourier_coeffs - Fourier Series coefficients in the format [base sine coeff, cosine1, sine1, cosine2, sine2,..., cosineN, sineN]

% upsampled_waveform - row vector containing the signal waveform (normalized to unit amplitude)
% time_vector [sec]
function [upsampled_waveform,time_vector,unnormalized_peak_amplitude] = construct_waveform(freq_i,fourier_coeffs)

var = num2cell(fourier_coeffs);
n_vars_minus_one = length(var)-1;
n_vars = size(var,2);

max_harmonic = ((n_vars_minus_one)/2)+1;
max_freq_component = freq_i*max_harmonic;

% determined empirically by plotting the following line with different N values: 
%       %LINE: figure; plot(sin(linspace(0,2*pi,N)),'.-')
% 13 points is clearly too coarse; 25 points is not bad; 50 point is awesome; 100 points is overkill
n_samples_per_max_freq_cycle = 50; 
period = 1/freq_i; % [sec]
dt = (1/(max_harmonic*freq_i))/n_samples_per_max_freq_cycle; % [sec]
time_vector = 0:dt:(period-dt); % [sec]

%%% handle first harmonic
harmonic_idx = 1;
flag_apply_sigma_correction = 0;
if (flag_apply_sigma_correction)
    correction_func = @(x) sinc(x/(max_harmonic+1));
else
    correction_func = @(x) 1;
end
signal = correction_func(harmonic_idx)*var{1,1}*sin(2*pi*freq_i*time_vector); % add up the first term

%%% handle higher harmonics
for n=2:2:(n_vars)
    harmonic_idx = (n/2)+1;
    signal = signal + correction_func(harmonic_idx)*(...
        var{1,n}*cos(2*pi*(harmonic_idx*freq_i)*time_vector) + var{1,n+1}*sin(2*pi*(harmonic_idx*freq_i)*time_vector));
end

upsampled_waveform = signal;
unnormalized_peak_amplitude = max(abs(upsampled_waveform));
upsampled_waveform = upsampled_waveform/unnormalized_peak_amplitude;
