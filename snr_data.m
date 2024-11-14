function snr = snr_data(no_signal, signal)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
snr = mean(std(signal))/mean(std(no_signal));
end