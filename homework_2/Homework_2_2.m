%% Homework 2 Part 2


close all, clc, clear all

%% Define colors

orange = [1 0.7 0];
orange2 = [0.5 0.6 0.5];
red = [0.8 0.1 0.1];

%% Load signal 1

[y,Fs] = audioread('music1.wav');
signal = y';
signal_data_points = length(y);
tr_piano=signal_data_points/Fs; % record time in seconds
time_vector = (1:length(y))/Fs;
signal_time = length(y)/Fs;

figure(1)

plot(time_vector,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (piano)');
p8 = audioplayer(y,Fs); playblocking(p8);

pause(0.000001)
%% Fourier of signal 1

frequencies_space = (2*pi/tr_piano)*[0:(signal_data_points/2 -1) -signal_data_points/2:-1];
frequencies_space_shifted = fftshift(frequencies_space);

transformed_signal = fft(signal);

figure(2)
plot(frequencies_space_shifted,fftshift(abs(transformed_signal)),  'Color', red);
xlabel('Frequncies');
ylabel('Amplitude');
title('Frequencies of Interest');
axis([-25000 25000 0 6000])

pause(0.000001)
%% Gauss filter
tau = signal_time/2;
a = 1;
gauss_filter = exp(-a*(time_vector-tau).^2);

figure(3)

plot(time_vector,gauss_filter, 'Linewidth', 2,  'Color', orange)
axis([0 signal_time 0 1.1])
xlabel('Time [sec]');
ylabel('Amplitude');
title('Gauss Filter');

pause(0.000001)
%% Animation signal 1

steps_per_second = 20;
steps = steps_per_second * tr_piano;
tau0 = 0;
a = 100;

figure(4)

for tau = 1:steps
    
    gauss_filter = exp(-a*(tr_piano - (tau0 + tau/steps_per_second)).^2);
    
    1*1
    
    vg = gauss_filter.*signal;
    
    1*2

    vgf = fft(vg);
    subplot(3,1,1)
    plot(frequencies_space_shifted,fftshift(abs(transformed_signal)),  'Color', red);
    axis([-10000 10000 0 6000])
    xlabel('Frequencies');
    ylabel('Amplitude');
    title('Frequencies of whole signal');
    
    1*3
    
    %subplot(3,1,2)
    %plot(time_vector,signal);
    hold on
    %plot(time_vector,gauss_filter, 'Linewidth', 2, 'Color', orange)
    %plot(time_vector,vg, 'Color', orange2)
    hold off
    xlabel('Time [sec]');
    ylabel('Amplitude');
    title('Gauss Filter');
    
    1*4
    
    %subplot(3,1,3)
    %plot(frequencies_space_shifted,fftshift(abs(vgf)),  'Color', red);
    %axis([-signal_data_points/2 signal_data_points/2 0 100])
    %xlabel('Frequencies');
    %ylabel('Amplitude');
    %title('Frequencies of filtered signal');
    
    
    
    pause(0.1)
    
end

pause(0.000001)
%% load signal 2

figure(5)
[y,Fs] = audioread('music2.wav');
tr_rec=length(y)/Fs; % record time in seconds
time_vector = (1:length(y))/Fs;
plot(time_vector,y);
xlabel('Time [sec]'); ylabel('Amplitude');
title('Mary had a little lamb (recorder)');
p8 = audioplayer(y,Fs); playblocking(p8);

