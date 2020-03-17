%% Homework 4

clear all, clc, close all

%% Otis McDonald (Hip Hop)

%% load test files

samples_per_song = 10;

number_of_songs = 5;

number_of_samples = samples_per_song * number_of_songs;

sample = cell(samples_per_song, 1);

for k = 1:samples_per_song
    
    sample{k} = [1500000+((k-1)*225000) 1500000+(k*225000)];
    
end

y_otis = cell(number_of_samples);
Fs_otis = cell(number_of_samples);

for k = 1:samples_per_song
    
    [y_otis{k}, Fs_otis{k}] = audioread("Behind_Closed_Doors(otis).mp3", sample{k});

    [y_otis{k+1*samples_per_song}, Fs_otis{k+1*samples_per_song}] = audioread("La_La_La(otis).mp3", sample{k});
    
    [y_otis{k+2*samples_per_song}, Fs_otis{k+2*samples_per_song}] = audioread("Mighty_Fine(otis).mp3", sample{k});
    
    [y_otis{k+3*samples_per_song}, Fs_otis{k+3*samples_per_song}] = audioread("Otis_McMusic(otis).mp3", sample{k});
    
    [y_otis{k+4*samples_per_song}, Fs_otis{k+4*samples_per_song}] = audioread("Sneaking_on_September(otis).mp3", sample{k});
    
end


%% transform data

signal_otis = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    signal_otis{k} = y_otis{k}(:, 1)';
    
end

%% define signal properties

data_points_otis = cell(number_of_samples, 1);
time_otis = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_points_otis{k} = length(signal_otis{k});
    time_otis{k} = data_points_otis{k}/Fs_otis{k};
    
end

%% define vectors

data_vector_otis = cell(number_of_samples, 1);
time_vector_otis = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_vector_otis{k} = 1:data_points_otis{k};
    time_vector_otis{k} = data_vector_otis{k}/Fs_otis{k};
    
end

%% sound plot

figure(1)

for k = 1:number_of_samples
    
    subplot(number_of_samples, 1, k)
    plot(time_vector_otis{k},signal_otis{k});
    axis off
    %axis([0 time_otis{k} -1 1])
    %xlabel('Time [sec]');
    %ylabel('Amplitude');
    %graph_name = sprintf('Signal of Interest, otis%d',k);
    %title(graph_name);
    
end

%% Fourier transform

fft_otis = cell(number_of_samples, 1);
frequencies_space_otis = cell(number_of_samples, 1);
frequencies_space_shifted_otis = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fft_otis{k} = fft(signal_otis{k});
    
    frequencies_space_otis{k} = (2*pi/time_otis{k}*[0:(data_points_otis{k}/2) -data_points_otis{k}/2:-1]);
    frequencies_space_shifted_otis{k} = fftshift(frequencies_space_otis{k});
    
    %for k = number_of_samples
        
    %    figure(9)

    %    plot(frequencies_space_shifted_otis{k},fftshift(abs(fft_otis{k})));
    %    xlabel('Frequncies');
    %    ylabel('Amplitude');
    %    graph_name = sprintf('Frequencies of Interest, otis%d',k);
    %    title(graph_name);
        %file_name = sprintf('otis_%d.png',k);
        %saveas(gcf, file_name)
        
    %end
    
    %figure(2)

    %plot(frequencies_space_shifted_otis{k},fftshift(abs(fft_otis{k})));
    %xlabel('Frequncies');
    %ylabel('Amplitude');
    %graph_name = sprintf('Frequencies of Interest, otis%d',k);
    %title(graph_name);
    %file_name = sprintf('otis_%d.png',k);
    %saveas(gcf, file_name)
    
end

%% preparing data

fourier_otis = cell(number_of_samples, 1);

spread_factor = 1000;

spread_variable = 5;

for k = 1:number_of_samples
    
    fourier_otis{k}(1, :) = frequencies_space_shifted_otis{k};
    fourier_otis{k}(2, :) = fftshift(abs(fft_otis{k}));
    
end

length_fourier_otis = cell(number_of_samples, 1);

sum_fourier_otis = cell(number_of_samples, 1);
 
max_fourier_otis = cell(number_of_samples, 1);
 
mean_fourier_otis = cell(number_of_samples, 1);

spread_fourier_otis = cell(number_of_samples, 1);

for j = 1 : number_of_samples
    
    length_fourier_otis{j} = length(fourier_otis{j});
    
    max_fourier_otis{j} = max(abs(fourier_otis{j}(1, j)));
    
    sum_fourier_otis{j} = 0;
    
    mean_fourier_otis{j} = 0;
    
    spread_fourier_otis{j} = 0;
    
    for k = 1 : length_fourier_otis{j}
        
        sum_fourier_otis{j} = sum_fourier_otis{j} + fourier_otis{j}(2, k);
        
        mean_fourier_otis{j} = mean_fourier_otis{j} + abs(fourier_otis{j}(1, k)) * fourier_otis{j}(2, k) / max_fourier_otis{j};
        
        if fourier_otis{j}(2, k) > (max_fourier_otis{j} / spread_factor)
            
            spread_fourier_otis{j} = spread_fourier_otis{j} + 1;
            
        end
        
    end
    
    sum_fourier_otis{j}  = sum_fourier_otis{j} / length_fourier_otis{j};
    
    mean_fourier_otis{j} = mean_fourier_otis{j} / length_fourier_otis{j};
    
    spread_fourier_otis{j} = (spread_fourier_otis{j})^(1/spread_variable);
    
end


%% E's Jammy Jams (Jazz & Blues)

%% load files

y_jammy = cell(number_of_samples);
Fs_jammy = cell(number_of_samples);

for k = 1:samples_per_song
    
    [y_jammy{k}, Fs_jammy{k}] = audioread("Happy_Birthday(jammy).mp3", sample{k});

    [y_jammy{k+1*samples_per_song}, Fs_jammy{k+1*samples_per_song}] = audioread("Maple_Leaf_Rag(jammy).mp3", sample{k});
    
    [y_jammy{k+2*samples_per_song}, Fs_jammy{k+2*samples_per_song}] = audioread("Minor_Blues_for_Booker(jammy).mp3", sample{k});
    
    [y_jammy{k+3*samples_per_song}, Fs_jammy{k+3*samples_per_song}] = audioread("Present_Day(jammy).mp3", sample{k});
    
    [y_jammy{k+4*samples_per_song}, Fs_jammy{k+4*samples_per_song}] = audioread("When_Johnny_Goes_Marching(jammy).mp3", sample{k});
    
end


%% transform data

signal_jammy = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    signal_jammy{k} = y_jammy{k}(:, 1)';
    
end

%% define signal properties

data_points_jammy = cell(number_of_samples, 1);
time_jammy = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_points_jammy{k} = length(signal_jammy{k});
    time_jammy{k} = data_points_jammy{k}/Fs_jammy{k};
    
end

%% define vectors

data_vector_jammy = cell(number_of_samples, 1);
time_vector_jammy = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_vector_jammy{k} = 1:data_points_jammy{k};
    time_vector_jammy{k} = data_vector_jammy{k}/Fs_jammy{k};
    
end

%% sound plot

figure(2)

for k = 1:number_of_samples
    
    subplot(number_of_samples, 1, k)
    plot(time_vector_jammy{k},signal_jammy{k});
    axis off
    %axis([0 time_jammy{k} -1 1])
    %xlabel('Time [sec]');
    %ylabel('Amplitude');
    %graph_name = sprintf('Signal of Interest, jammy%d',k);
    %title(graph_name);
    
end

%% Fourier transform

fft_jammy = cell(number_of_samples, 1);
frequencies_space_jammy = cell(number_of_samples, 1);
frequencies_space_shifted_jammy = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fft_jammy{k} = fft(signal_jammy{k});
    
    frequencies_space_jammy{k} = (2*pi/time_jammy{k}*[0:(data_points_jammy{k}/2) -data_points_jammy{k}/2:-1]);
    frequencies_space_shifted_jammy{k} = fftshift(frequencies_space_jammy{k});
    
    %for k = number_of_samples
        
    %    figure(10)

    %    plot(frequencies_space_shifted_jammy{k},fftshift(abs(fft_jammy{k})));
    %    xlabel('Frequncies');
    %    ylabel('Amplitude');
    %    graph_name = sprintf('Frequencies of Interest, jammy%d',k);
    %    title(graph_name);
        %file_name = sprintf('jammy_%d.png',k);
        %saveas(gcf, file_name)
        
    %end
    
end

%% preparing data

fourier_jammy = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fourier_jammy{k}(1, :) = frequencies_space_shifted_jammy{k};
    fourier_jammy{k}(2, :) = fftshift(abs(fft_jammy{k}));
    
end

length_fourier_jammy = cell(number_of_samples, 1);

sum_fourier_jammy = cell(number_of_samples, 1);
 
max_fourier_jammy = cell(number_of_samples, 1);
 
mean_fourier_jammy = cell(number_of_samples, 1);

spread_fourier_jammy = cell(number_of_samples, 1);

for j = 1 : number_of_samples
    
    length_fourier_jammy{j} = length(fourier_jammy{j});
    
    max_fourier_jammy{j} = max(abs(fourier_jammy{j}(1, j)));
    
    sum_fourier_jammy{j} = 0;
    
    mean_fourier_jammy{j} = 0;
    
    spread_fourier_jammy{j} = 0;
    
    for k = 1 : length_fourier_jammy{j}
        
        sum_fourier_jammy{j} = sum_fourier_jammy{j} + fourier_jammy{j}(2, k);
        
        mean_fourier_jammy{j} = mean_fourier_jammy{j} + abs(fourier_jammy{j}(1, k)) * fourier_jammy{j}(2, k) / max_fourier_jammy{j};
        
        if fourier_jammy{j}(2, k) > (max_fourier_jammy{j} / spread_factor)
            
            spread_fourier_jammy{j} = spread_fourier_jammy{j} + 1;
            
        end
        
    end
    
    sum_fourier_jammy{j}  = sum_fourier_jammy{j} / length_fourier_jammy{j};
    
    mean_fourier_jammy{j} = mean_fourier_jammy{j} / length_fourier_jammy{j};
    
    spread_fourier_jammy{j} = (spread_fourier_jammy{j})^(1/spread_variable);
    
end


%% Max MacFerren (Classic)

%% load files

y_max = cell(number_of_samples);
Fs_max = cell(number_of_samples);

for k = 1:samples_per_song
    
    [y_max{k}, Fs_max{k}] = audioread( "All_I_Can_Do_Is_This(max).mp3", sample{k});

    [y_max{k+1*samples_per_song}, Fs_max{k+1*samples_per_song}] = audioread("Come_With_Some_Funk(max).mp3", sample{k});
    
    [y_max{k+2*samples_per_song}, Fs_max{k+2*samples_per_song}] = audioread("Mom_s_House(max).mp3", sample{k});
    
    [y_max{k+3*samples_per_song}, Fs_max{k+3*samples_per_song}] = audioread("Sense_of_Humor(max).mp3", sample{k});
    
    [y_max{k+4*samples_per_song}, Fs_max{k+4*samples_per_song}] = audioread("Waterfront_Property(max).mp3", sample{k});
    
end


%% transform data

signal_max = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    signal_max{k} = y_max{k}(:, 1)';
    
end

%% define signal properties

data_points_max = cell(number_of_samples, 1);
time_max = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_points_max{k} = length(signal_max{k});
    time_max{k} = data_points_max{k}/Fs_max{k};
    
end

%% define vectors

data_vector_max = cell(number_of_samples, 1);
time_vector_max = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_vector_max{k} = 1:data_points_max{k};
    time_vector_max{k} = data_vector_max{k}/Fs_max{k};
    
end

%% sound plot

figure(3)

for k = 1:number_of_samples
    
    subplot(number_of_samples, 1, k)
    plot(time_vector_max{k},signal_max{k});
    axis off
    %axis([0 time_max{k} -1 1])
    %xlabel('Time [sec]');
    %ylabel('Amplitude');
    %graph_name = sprintf('Signal of Interest, max%d',k);
    %title(graph_name);
    
end

%% Fourier transform

fft_max = cell(number_of_samples, 1);
frequencies_space_max = cell(number_of_samples, 1);
frequencies_space_shifted_max = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fft_max{k} = fft(signal_max{k});
    
    frequencies_space_max{k} = (2*pi/time_max{k}*[0:(data_points_max{k}/2) -data_points_max{k}/2:-1]);
    frequencies_space_shifted_max{k} = fftshift(frequencies_space_max{k});
    
    %for k = number_of_samples
        
    %    figure(11)

    %    plot(frequencies_space_shifted_max{k},fftshift(abs(fft_max{k})));
    %    xlabel('Frequncies');
    %    ylabel('Amplitude');
    %    graph_name = sprintf('Frequencies of Interest, max%d',k);
    %    title(graph_name);
        %file_name = sprintf('max_%d.png',k);
        %saveas(gcf, file_name)
        
    %end
    
end

%% preparing data

fourier_max = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fourier_max{k}(1, :) = frequencies_space_shifted_max{k};
    fourier_max{k}(2, :) = fftshift(abs(fft_max{k}));
    
end

length_fourier_max = cell(number_of_samples, 1);

sum_fourier_max = cell(number_of_samples, 1);
 
max_fourier_max = cell(number_of_samples, 1);
 
mean_fourier_max = cell(number_of_samples, 1);

spread_fourier_max = cell(number_of_samples, 1);

for j = 1 : number_of_samples
    
    length_fourier_max{j} = length(fourier_max{j});
    
    max_fourier_max{j} = max(abs(fourier_max{j}(1, j)));
    
    sum_fourier_max{j} = 0;
    
    mean_fourier_max{j} = 0;
    
    spread_fourier_max{j} = 0;
    
    for k = 1 : length_fourier_max{j}
        
        sum_fourier_max{j} = sum_fourier_max{j} + fourier_max{j}(2, k);
        
        mean_fourier_max{j} = mean_fourier_max{j} + abs(fourier_max{j}(1, k)) * fourier_max{j}(2, k) / max_fourier_max{j};
        
        if fourier_max{j}(2, k) > (max_fourier_max{j} / spread_factor)
            
            spread_fourier_max{j} = spread_fourier_max{j} + 1;
            
        end
        
    end
    
    sum_fourier_max{j}  = sum_fourier_max{j} / length_fourier_max{j};
    
    mean_fourier_max{j} = mean_fourier_max{j} / length_fourier_max{j};
    
    spread_fourier_max{j} = (spread_fourier_max{j})^(1/spread_variable);
    
end


%% All together

%% plot data

figure(4)

sum_fourier_otis_plot = zeros(1, number_of_samples);

mean_fourier_otis_plot = zeros(1, number_of_samples);

spread_fourier_otis_plot = zeros(1, number_of_samples);

for k = 1:number_of_samples
    
    sum_fourier_otis_plot(k) = sum_fourier_otis{k};
    
    mean_fourier_otis_plot(k) = mean_fourier_otis{k};
    
    spread_fourier_otis_plot(k) = spread_fourier_otis{k};
    
end

%plot(sum_fourier_otis_plot, mean_fourier_otis_plot, 'bo')
scatter3(sum_fourier_otis_plot, mean_fourier_otis_plot, spread_fourier_otis_plot, 'blue')
hold on
axis vis3d


sum_fourier_jammy_plot = zeros(1, number_of_samples);

mean_fourier_jammy_plot = zeros(1, number_of_samples);

spread_fourier_jammy_plot = zeros(1, number_of_samples);

for k = 1:number_of_samples
    
    sum_fourier_jammy_plot(k) = sum_fourier_jammy{k};
    
    mean_fourier_jammy_plot(k) = mean_fourier_jammy{k};
    
    spread_fourier_jammy_plot(k) = spread_fourier_jammy{k};
    
end

%plot(sum_fourier_jammy_plot, mean_fourier_jammy_plot, 'ro')
scatter3(sum_fourier_jammy_plot, mean_fourier_jammy_plot, spread_fourier_jammy_plot, 'red')


sum_fourier_max_plot = zeros(1, number_of_samples);

mean_fourier_max_plot = zeros(1, number_of_samples);

spread_fourier_max_plot = zeros(1, number_of_samples);

for k = 1:number_of_samples
    
    sum_fourier_max_plot(k) = sum_fourier_max{k};
    
    mean_fourier_max_plot(k) = mean_fourier_max{k};
    
    spread_fourier_max_plot(k) = spread_fourier_max{k};
    
end

%plot(sum_fourier_max_plot, mean_fourier_max_plot, 'ko')
scatter3(sum_fourier_max_plot, mean_fourier_max_plot, spread_fourier_max_plot, 'black')

%% mean plot

mean_otis(1, 1) = mean(sum_fourier_otis_plot);

mean_otis(1, 2) = mean(mean_fourier_otis_plot);

mean_otis(1, 3) = mean(spread_fourier_otis_plot);

scatter3(mean_otis(1, 1), mean_otis(1, 2), mean_otis(1, 3), 'yellow')
%plot(mean_otis(1, 1), mean_otis(1, 2), 'y*')


mean_jammy(1, 1) = mean(sum_fourier_jammy_plot);

mean_jammy(1, 2) = mean(mean_fourier_jammy_plot);

mean_jammy(1, 3) = mean(spread_fourier_jammy_plot);

scatter3(mean_jammy(1, 1), mean_jammy(1, 2), mean_jammy(1, 3), 'yellow')
%plot(mean_jammy(1, 1), mean_jammy(1, 2), 'y*')


mean_max(1, 1) = mean(sum_fourier_max_plot);

mean_max(1, 2) = mean(mean_fourier_max_plot);

mean_max(1, 3) = mean(spread_fourier_max_plot);

scatter3(mean_max(1, 1), mean_max(1, 2), mean_max(1, 3), 'yellow')
%plot(mean_max(1, 1), mean_max(1, 2), 'y*')


%% covariance

figure(5)

zero_mean_otis(1, :) = sum_fourier_otis_plot - mean_otis(1, 1);

zero_mean_otis(2, :) = mean_fourier_otis_plot - mean_otis(1, 2);

zero_mean_otis(3, :) = spread_fourier_otis_plot - mean_otis(1, 3);

scatter3(zero_mean_otis(1, :), zero_mean_otis(2, :), zero_mean_otis(3, :), 'blue')
axis vis3d

%plot(zero_mean_otis(1, :), zero_mean_otis(2, :), 'bo')

covarience_sum_otis = cov(zero_mean_otis(1, :));

covariance_mean_otis = cov(zero_mean_otis(2, :));

covariance_spread_otis = cov(zero_mean_otis(3, :));


figure(6)

zero_mean_jammy(1, :) = sum_fourier_jammy_plot - mean_jammy(1, 1);

zero_mean_jammy(2, :) = mean_fourier_jammy_plot - mean_jammy(1, 2);

zero_mean_jammy(3, :) = spread_fourier_jammy_plot - mean_jammy(1, 3);

scatter3(zero_mean_jammy(1, :), zero_mean_jammy(2, :), zero_mean_jammy(3, :), 'red')
axis vis3d

%plot(zero_mean_jammy(1, :), zero_mean_jammy(2, :), 'ro')

covarience_sum_jammy = cov(zero_mean_jammy(1, :));

covariance_mean_jammy = cov(zero_mean_jammy(2, :));

covariance_spread_jammy = cov(zero_mean_jammy(3, :));


figure(7)

zero_mean_max(1, :) = sum_fourier_max_plot - mean_max(1, 1);

zero_mean_max(2, :) = mean_fourier_max_plot - mean_max(1, 2);

zero_mean_max(3, :) = spread_fourier_max_plot - mean_max(1, 3);

scatter3(zero_mean_max(1, :), zero_mean_max(2, :), zero_mean_max(3, :), 'black')
axis vis3d

%plot(zero_mean_max(1, :), zero_mean_max(2, :), 'ko')

covarience_sum_max = cov(zero_mean_max(1, :));

covariance_mean_max = cov(zero_mean_max(2, :));

covariance_spread_max = cov(zero_mean_max(3, :));

%% SVD Plot

figure(8)

[U_otis, S_otis, V_otis] = svd(zero_mean_otis);

[m_otis, n_otis]=size(zero_mean_otis);

S2_otis=S_otis(1:m_otis,1:m_otis); Values_otis=S2_otis^2/n_otis;  %Eigenvalues

Vectors_otis=U_otis;                     %Eigenvectors

Amplitudes_otis=S_otis*conj(V_otis');         %Orthognal amplitudes

scatter_otis = [zero_mean_otis(1,:) + mean_otis(1, 1);
                zero_mean_otis(2,:) + mean_otis(1, 2);
                zero_mean_otis(3, :) + mean_otis(1, 3)];

scatter3(scatter_otis(1, :), scatter_otis(2, :), scatter_otis(3, :),'b.');
hold on
axis vis3d
title('Train Data');

line1_otis = [(Values_otis(1,1)*[0 ; Vectors_otis(1,1)] + mean_otis(1, 1))';
              (Values_otis(1,1)*[0 ; Vectors_otis(2,1)] + mean_otis(1, 2))';
              (Values_otis(1,1)*[0 ; Vectors_otis(3,1)] + mean_otis(1, 3))'];
          
          
r1_otis = [line1_otis(1, 2) - line1_otis(1, 1), line1_otis(2, 2) - line1_otis(2, 1), line1_otis(3, 2) - line1_otis(3, 1)];

%r2_otis = [line2_otis(1, 2) - line2_otis(1, 1), line2_otis(2, 2) - line2_otis(2, 1), line2_otis(3, 2) - line2_otis(3, 1)];

r_otis = abs(r1_otis);

p_otis(:, 1) = mean_otis' - 0.5 * r_otis';

p_otis(:, 2) = mean_otis' + 0.5 * r_otis';

plot3(p_otis(1, :), p_otis(2, :), p_otis(3, :), 'b')


[U_jammy, S_jammy, V_jammy] = svd(zero_mean_jammy);

[m_jammy, n_jammy]=size(zero_mean_jammy);

S2_jammy=S_jammy(1:m_jammy,1:m_jammy); Values_jammy=S2_jammy^2/n_jammy;  %Eigenvalues

Vectors_jammy=U_jammy;                     %Eigenvectors

Amplitudes_jammy=S_jammy*conj(V_jammy');         %Orthognal amplitudes

scatter_jammy = [zero_mean_jammy(1,:) + mean_jammy(1, 1);
                zero_mean_jammy(2,:) + mean_jammy(1, 2);
                zero_mean_jammy(3, :) + mean_jammy(1, 3)];
            
scatter3(scatter_jammy(1, :), scatter_jammy(2, :), scatter_jammy(3, :),'r.');

line1_jammy = [(Values_jammy(1,1)*[0 ; Vectors_jammy(1,1)] + mean_jammy(1, 1))';
               (Values_jammy(1,1)*[0 ; Vectors_jammy(2,1)] + mean_jammy(1, 2))';
               (Values_jammy(1,1)*[0 ; Vectors_jammy(3,1)] + mean_jammy(1, 3))'];
           
r1_jammy = [line1_jammy(1, 2) - line1_jammy(1, 1), line1_jammy(2, 2) - line1_jammy(2, 1), line1_jammy(3, 2) - line1_jammy(3, 1)];

%r2_jammy = [line2_jammy(1, 2) - line2_jammy(1, 1), line2_jammy(2, 2) - line2_jammy(2, 1), line2_jammy(3, 2) - line2_jammy(3, 1)];

r_jammy = abs(r1_jammy);

p_jammy(:, 1) = mean_jammy' - r_jammy';

p_jammy(:, 2) = mean_jammy' + r_jammy';

plot3(p_jammy(1, :), p_jammy(2, :), p_jammy(3, :), 'r')

[U_max, S_max, V_max] = svd(zero_mean_max);

[m_max, n_max]=size(zero_mean_max);

S2_max=S_max(1:m_max,1:m_max); Values_max=S2_max^2/n_max;  %Eigenvalues

Vectors_max=U_max;                     %Eigenvectors

Amplitudes_max=S_max*conj(V_max');         %Orthognal amplitudes

scatter_max = [zero_mean_max(1,:) + mean_max(1, 1);
                zero_mean_max(2,:) + mean_max(1, 2);
                zero_mean_max(3, :) + mean_max(1, 3)];
            
scatter3(scatter_max(1, :), scatter_max(2, :), scatter_max(3, :),'k.');

line1_max = [(Values_max(1,1)*[0 ; Vectors_max(1,1)] + mean_max(1, 1))';
             (Values_max(1,1)*[0 ; Vectors_max(2,1)] + mean_max(1, 2))';
             (Values_max(1,1)*[0 ; Vectors_max(3,1)] + mean_max(1, 3))'];
         
line2_max = [(Values_max(1,1)*[0 ; Vectors_max(1,2)] + mean_max(1, 1))';
             (Values_max(1,1)*[0 ; Vectors_max(2,2)] + mean_max(1, 2))';
             (Values_max(1,1)*[0 ; Vectors_max(3,2)] + mean_max(1, 3))'];

r1_max = [line1_max(1, 2) - line1_max(1, 1), line1_max(2, 2) - line1_max(2, 1), line1_max(3, 2) - line1_max(3, 1)];

r2_max = [line2_max(1, 2) - line2_max(1, 1), line2_max(2, 2) - line2_max(2, 1), line2_max(3, 2) - line2_max(3, 1)];

r_max = abs(r1_max);

p_max(:, 1) = mean_max' - 1.5 * r_max';

p_max(:, 2) = mean_max' + r_max';

plot3(p_max(1, :), p_max(2, :), p_max(3, :), 'k')

%% Chosing sector

weight1 = 10;

weight2 = 1;


d_otis_right = point_to_line_distance(scatter_otis', p_otis(:, 1)', p_otis(:, 2)');

for k = 1: length(d_otis_right)
    
    d_otis_right(k, :) = weight1 * d_otis_right(k, :) + weight2 * norm(abs(scatter_otis(:, k)' - mean_otis));
     
end

d_otis_wrong1 = point_to_line_distance(scatter_jammy', p_otis(:, 1)', p_otis(:, 2)');

d_otis_wrong2 = point_to_line_distance(scatter_max', p_otis(:, 1)', p_otis(:, 2)');

d_otis_wrong = zeros(length(d_otis_wrong1) + length(d_otis_wrong2), 1);

for k = 1 : length(d_otis_wrong1)
    
    d_otis_wrong(k, 1) = weight1 * d_otis_wrong1(k) + weight2 * norm(abs(scatter_jammy(:, k)' - mean_otis));
    
end

for k = 1 : k
    
    d_otis_wrong(k + length(d_otis_wrong1), 1) = weight1 * d_otis_wrong2(k) + weight2 * norm(abs(scatter_max(:, k)' - mean_otis));
    
end

sector_otis = max(d_otis_right);

next_sector_otis = min(d_otis_wrong);


d_jammy_right = point_to_line_distance(scatter_jammy', p_jammy(:, 1)', p_jammy(:, 2)');

for k = 1: length(d_jammy_right)
    
    d_jammy_right(k, :) = weight1 * d_jammy_right(k, :) + weight2 * norm(abs(scatter_jammy(:, k)' - mean_jammy));
    
end

d_jammy_wrong1 = point_to_line_distance(scatter_otis', p_jammy(:, 1)', p_jammy(:, 2)');

d_jammy_wrong2 = point_to_line_distance(scatter_max', p_jammy(:, 1)', p_jammy(:, 2)');

d_jammy_wrong = zeros(length(d_jammy_wrong1) + length(d_jammy_wrong2), 1);

for k = 1 : length(d_jammy_wrong1)
    
    d_jammy_wrong(k, 1) = weight1 * d_jammy_wrong1(k)  + weight2 * norm(abs(scatter_otis(:, k)' - mean_jammy));
    
end

for k = 1 : k
    
    d_jammy_wrong(k + length(d_jammy_wrong1), 1) = weight1 * d_jammy_wrong2(k) + weight2 * norm(abs(scatter_max(:, k)' - mean_jammy));
    
end

sector_jammy = max(d_jammy_right);

next_sector_jammy = min(d_jammy_wrong);


d_max_right = point_to_line_distance(scatter_max', p_max(:, 1)', p_max(:, 2)');

for k = 1 : length(d_max_right)
    
    d_max_right(k, :) = weight1 * d_max_right(k, :) + weight2 * norm(abs(scatter_max(:, k)' - mean_max));
    
end

d_max_wrong1 = point_to_line_distance(scatter_otis', p_max(:, 1)', p_max(:, 2)');

d_max_wrong2 = point_to_line_distance(scatter_jammy', p_max(:, 1)', p_max(:, 2)');

d_max_wrong = zeros(length(d_max_wrong1) + length(d_max_wrong2), 1);

for k = 1 : length(d_max_wrong1)
    
    d_max_wrong(k, 1) = weight1 * d_max_wrong1(k) + weight2 * norm(abs(scatter_otis(:, k)' - mean_max));
    
end

for k = 1 : k
    
    d_max_wrong(k + length(d_max_wrong1), 1) = weight1 * d_max_wrong2(k) + weight2 * norm(abs(scatter_jammy(:, k)' - mean_max));
    
end

sector_max = max(d_max_right);

next_sector_max = min(d_max_wrong);

% Building borders

next_sector = [next_sector_otis, next_sector_jammy, next_sector_max];

border_otis = (sector_otis + next_sector_otis);

border_jammy = (sector_jammy + next_sector_jammy) / 2;

border_max = (sector_max + next_sector_max);

%% Test with new data

%% load test files

samples_per_song = 10;

number_of_songs = 6;

number_of_samples_per_genre = samples_per_song * 2;

number_of_samples = samples_per_song * number_of_songs;

sample = cell(samples_per_song, 1);

for k = 1:samples_per_song
    
    sample{k} = [2000000+((k-1)*225000) 2000000+(k*225000)];
    
end

y_test = cell(number_of_samples);
Fs_test = cell(number_of_samples);

for k = 1:samples_per_song
    
    [y_test{k}, Fs_test{k}] = audioread("Complicate_ya(hiphop).mp3", sample{k});

    [y_test{k+1*samples_per_song}, Fs_test{k+1*samples_per_song}] = audioread("Here_If_You_re_Going(hiphop).mp3", sample{k});
    
    [y_test{k+2*samples_per_song}, Fs_test{k+2*samples_per_song}] = audioread("Present_Day(jazz).mp3", sample{k});
    
    [y_test{k+3*samples_per_song}, Fs_test{k+3*samples_per_song}] = audioread("Tango_Bango(jazz).mp3", sample{k});
    
    [y_test{k+4*samples_per_song}, Fs_test{k+4*samples_per_song}] = audioread("Cosmic_Love(Electro).mp3", sample{k});
    
    [y_test{k+5*samples_per_song}, Fs_test{k+5*samples_per_song}] = audioread("Tired_Break(Electro).mp3", sample{k});
    
end

%% transform data

signal_test = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    signal_test{k} = y_test{k}(:, 1)';
    
end

%% define signal properties

data_points_test = cell(number_of_samples, 1);
time_test = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_points_test{k} = length(signal_test{k});
    time_test{k} = data_points_test{k}/Fs_test{k};
    
end

%% define vectors

data_vector_test = cell(number_of_samples, 1);
time_vector_test = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    data_vector_test{k} = 1:data_points_test{k};
    time_vector_test{k} = data_vector_test{k}/Fs_test{k};
    
end

%% Fourier transform

fft_test = cell(number_of_samples, 1);
frequencies_space_test = cell(number_of_samples, 1);
frequencies_space_shifted_test = cell(number_of_samples, 1);

for k = 1:number_of_samples
    
    fft_test{k} = fft(signal_test{k});
    
    frequencies_space_test{k} = (2*pi/time_test{k}*[0:(data_points_test{k}/2) -data_points_test{k}/2:-1]);
    frequencies_space_shifted_test{k} = fftshift(frequencies_space_test{k});
    
end

%% preparing data

fourier_test = cell(number_of_samples, 1);

spread_factor = 1000;

spread_variable = 5;

for k = 1:number_of_samples
    
    fourier_test{k}(1, :) = frequencies_space_shifted_test{k};
    fourier_test{k}(2, :) = fftshift(abs(fft_test{k}));
    
end

length_fourier_test = cell(number_of_samples, 1);

sum_fourier_test = cell(number_of_samples, 1);
 
max_fourier_test = cell(number_of_samples, 1);
 
mean_fourier_test = cell(number_of_samples, 1);

spread_fourier_test = cell(number_of_samples, 1);

for j = 1 : number_of_samples
    
    length_fourier_test{j} = length(fourier_test{j});
    
    max_fourier_test{j} = max(abs(fourier_test{j}(1, j)));
    
    sum_fourier_test{j} = 0;
    
    mean_fourier_test{j} = 0;
    
    spread_fourier_test{j} = 0;
    
    for k = 1 : length_fourier_test{j}
        
        sum_fourier_test{j} = sum_fourier_test{j} + fourier_test{j}(2, k);
        
        mean_fourier_test{j} = mean_fourier_test{j} + abs(fourier_test{j}(1, k)) * fourier_test{j}(2, k) / max_fourier_test{j};
        
        if fourier_test{j}(2, k) > (max_fourier_test{j} / spread_factor)
            
            spread_fourier_test{j} = spread_fourier_test{j} + 1;
            
        end
        
    end
    
    sum_fourier_test{j}  = sum_fourier_test{j} / length_fourier_test{j};
    
    mean_fourier_test{j} = mean_fourier_test{j} / length_fourier_test{j};
    
    spread_fourier_test{j} = (spread_fourier_test{j})^(1/spread_variable);
    
end

%% plot data

figure(9)

plot3(p_otis(1, :), p_otis(2, :), p_otis(3, :), 'b')
hold on
title('Test Data')
plot3(p_jammy(1, :), p_jammy(2, :), p_jammy(3, :), 'r')
plot3(p_max(1, :), p_max(2, :), p_max(3, :), 'k')

sum_fourier_test_plot = zeros(1, number_of_samples);

mean_fourier_test_plot = zeros(1, number_of_samples);

spread_fourier_test_plot = zeros(1, number_of_samples);

for k = 1:number_of_samples
    
    sum_fourier_test_plot(k) = sum_fourier_test{k};
    
    mean_fourier_test_plot(k) = mean_fourier_test{k};
    
    spread_fourier_test_plot(k) = spread_fourier_test{k};
    
end

for k = 1 : number_of_samples_per_genre
    
    scatter3(sum_fourier_test_plot(k), mean_fourier_test_plot(k), spread_fourier_test_plot(k), 'b.')
    
end

for k = 1 : number_of_samples_per_genre
    
    scatter3(sum_fourier_test_plot(k + number_of_samples_per_genre), mean_fourier_test_plot(k + number_of_samples_per_genre), spread_fourier_test_plot(k + number_of_samples_per_genre), 'r.')
    
end

for k = 1 : number_of_samples_per_genre
    
    scatter3(sum_fourier_test_plot(k + 2*number_of_samples_per_genre), mean_fourier_test_plot(k + 2*number_of_samples_per_genre), spread_fourier_test_plot(k + 2*number_of_samples_per_genre), 'k.')
    
end

axis vis3d

%% Classify

scatter_hiphop = zeros(3, number_of_samples_per_genre);

scatter_jazz = zeros(3, number_of_samples_per_genre);

scatter_electro = zeros(3, number_of_samples_per_genre);

for k = 1 : number_of_samples_per_genre
    
    scatter_hiphop(:, k) =   [sum_fourier_test_plot(k);
                        mean_fourier_test_plot(k);
                        spread_fourier_test_plot(k)];
    
end


for k = 1 : number_of_samples_per_genre
    
    scatter_jazz(:, k) =   [sum_fourier_test_plot(k + number_of_samples_per_genre);
                        mean_fourier_test_plot(k + number_of_samples_per_genre);
                        spread_fourier_test_plot(k + number_of_samples_per_genre)];
    
end

for k = 1 : number_of_samples_per_genre
    
    scatter_electro(:, k) =   [sum_fourier_test_plot(k + 2 * number_of_samples_per_genre);
                        mean_fourier_test_plot(k + 2 * number_of_samples_per_genre);
                        spread_fourier_test_plot(k + 2 * number_of_samples_per_genre)];
    
end

hiphop_right = 0;

hiphop_wrong = 0;

not_hiphop_right = 0;

not_hiphop_wrong = 0;

for k = 1 : number_of_samples_per_genre
    
    d_otis_right(k, :) = weight1 * point_to_line_distance(scatter_hiphop(:, k)', p_otis(:, 1)', p_otis(:, 2)') + weight2 * norm(abs(scatter_hiphop(:, k)' - mean_otis));
    
    if d_otis_right(k, :) <= border_otis
        
        hiphop_right = hiphop_right + 1;
        
    else
        
        hiphop_wrong = hiphop_wrong + 1;
        
    end
    
    d_otis_wrong1(k, :) = weight1 * point_to_line_distance(scatter_jazz(:, k)', p_otis(:, 1)', p_otis(:, 2)') + weight2 * norm(abs(scatter_jazz(:, k)' - mean_otis));
    
    if d_otis_wrong1(k, :) > border_otis
        
        not_hiphop_right = not_hiphop_right + 1;
        
    else
        
        not_hiphop_wrong = not_hiphop_wrong + 1;
        
    end
    
    d_otis_wrong2(k, :) = weight1 * point_to_line_distance(scatter_electro(:, k)', p_otis(:, 1)', p_otis(:, 2)') + weight2 * norm(abs(scatter_electro(:, k)' - mean_otis));
    
    if d_otis_wrong2(k, :) > border_otis
        
        not_hiphop_right = not_hiphop_right + 1;
        
    else
        
        not_hiphop_wrong = not_hiphop_wrong + 1;
        
    end
    
end


jazz_right = 0;

jazz_wrong = 0;

not_jazz_right = 0;

not_jazz_wrong = 0;

for k = 1 : number_of_samples_per_genre
    
    d_jammy_right(k, :) = weight1 * point_to_line_distance(scatter_jazz(:, k)', p_jammy(:, 1)', p_jammy(:, 2)') + weight2 * norm(abs(scatter_jazz(:, k)' - mean_jammy));
    
    if d_jammy_right(k, :) <= border_jammy
        
        jazz_right = jazz_right + 1;
        
    else
        
        jazz_wrong = jazz_wrong + 1;
        
    end
    
    d_jammy_wrong1(k, :) = weight1 * point_to_line_distance(scatter_hiphop(:, k)', p_jammy(:, 1)', p_jammy(:, 2)') + weight2 * norm(abs(scatter_hiphop(:, k)' - mean_jammy));
    
    if d_jammy_wrong1(k, :) > border_jammy
        
        not_jazz_right = not_jazz_right + 1;
        
    else
        
        not_jazz_wrong = not_jazz_wrong + 1;
        
    end
    
    d_jammy_wrong2(k, :) = weight1 * point_to_line_distance(scatter_electro(:, k)', p_jammy(:, 1)', p_jammy(:, 2)') + weight2 * norm(abs(scatter_electro(:, k)' - mean_jammy));
    
    if d_jammy_wrong2(k, :) > border_jammy
        
        not_jazz_right = not_jazz_right + 1;
        
    else
        
        not_jazz_wrong = not_jazz_wrong + 1;
        
    end
    
end




electro_right = 0;

electro_wrong = 0;

not_electro_right = 0;

not_electro_wrong = 0;

for k = 1 : number_of_samples_per_genre
    
    d_max_right(k, :) = weight1 * point_to_line_distance(scatter_electro(:, k)', p_max(:, 1)', p_max(:, 2)') + weight2 * norm(abs(scatter_electro(:, k)' - mean_max));
    
    if d_max_right(k, :) <= border_max
        
        electro_right = electro_right + 1;
        
    else
        
        electro_wrong = electro_wrong + 1;
        
    end
    
    d_max_wrong1(k, :) = weight1 * point_to_line_distance(scatter_hiphop(:, k)', p_max(:, 1)', p_max(:, 2)') + weight2 * norm(abs(scatter_hiphop(:, k)' - mean_max));
    
    if d_max_wrong1(k, :) > border_max
        
        not_electro_right = not_electro_right + 1;
        
    else
        
        not_electro_wrong = not_electro_wrong + 1;
        
    end
    
    d_max_wrong2(k, :) = weight1 * point_to_line_distance(scatter_jazz(:, k)', p_max(:, 1)', p_max(:, 2)') + weight2 * norm(abs(scatter_jazz(:, k)' - mean_max));
    
    if d_max_wrong2(k, :) > border_max
        
        not_electro_right = not_electro_right + 1;
        
    else
        
        not_electro_wrong = not_electro_wrong + 1;
        
    end
    
end

quote_hiphop_right = hiphop_right / number_of_samples_per_genre

quote_not_hiphop_right = not_hiphop_right / (2 * number_of_samples_per_genre)


quote_jazz_right = jazz_right / number_of_samples_per_genre

quote_not_jazz_right = not_jazz_right / (2 * number_of_samples_per_genre)


quote_electro_right = electro_right / number_of_samples_per_genre

quote_not_electro_right = not_electro_right / (2 * number_of_samples_per_genre)



%% Functions

function distance=point_to_line_distance(pt, v1, v2)

if nargin~=3
    error('HJW:point_to_line_distance:nargin',...
        'Incorrect number of inputs, expected 3.');
end
if ~isnumeric(pt) || ~any(size(pt,2)==[2 3]) || any(size(pt)==0)
    error('HJW:point_to_line_distance:pt_type_size',...
        'First input (pt) is not numeric or has an incorrect shape.')
end
if ~isnumeric(v1) || numel(v1)~=size(pt,2)
    error('HJW:point_to_line_distance:v_type_size',...
        ['Second input (v1) is not numeric or has an incorrect ',...
        'size.' char(10) 'Expected 1x3 or 1x2, which should match ',...
        'the first input.']) %#ok<CHARTEN>
end
if ~isnumeric(v2) || numel(v2)~=size(pt,2)
    error('HJW:point_to_line_distance:v_type_size',['Third input (v2) ',...
        'is not numeric or has an incorrect size.' char(10) 'Expected ',...
        '1x3 or 1x2, which should match the first input.']) %#ok<CHARTEN>
end
v1=v1(:)';
if length(v1)==2,v1(3)=0;end
v2=v2(:)';
if length(v2)==2,v2(3)=0;end
if size(pt,2)==2,pt(1,3)=0;end
v1_ = repmat(v1,size(pt,1),1);
v2_ = repmat(v2,size(pt,1),1);

a = v1_ - v2_;
b = pt - v2_;
distance = sqrt(sum(cross(a,b,2).^2,2)) ./ sqrt(sum(a.^2,2));


end