clear variables; close all; clc

%Setting the path to the resources depending on the current OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-bsp-2023/src/Data/';
else
    path = 'Data\';
end

%Loading the data files into a dir struct
myFiles = dir(strcat(path, '*.mat'));

%Loading the file containing the specifics about the EEG
%TODO: clean
tmp = load(strcat(path, 'chanlocs.mat'));
chanlocs = tmp.chanlocs;

%Iterating through the activity data (skipping the first element)
for tmp = (length(myFiles) - 1)/2:-1:1
    filePath = strcat(path, myFiles(tmp * 2 - 1).name); %Odd spacing w/ given data
    myWorkData(tmp) = load(filePath);
end

%Iterating through the rest data
for tmp = (length(myFiles) - 1)/2:-1:1
    filePath = strcat(path, myFiles(tmp * 2).name); %Even spacing w/ given data
    myRestData(tmp) = load(filePath);
end

%Expliciting the struct attributes in order to iterate through them
%TODO: see if chanlocs order is the same as the signal one, use chanlocs in
%case it is
electrodes = ["C3","C4","CZ","F3","F4","F7","F8","FP1","FP2","FZ","O1","O2","P3","P4","PZ","T3","T4","T5","T6"];

Fs = 500;
% t = (1:30*Fs)/Fs;

%Defining the frequency scale
hz = linspace(0, Fs/2, floor(Fs/2) + 1);

%Expliciting the number of samples to ignore
sigErr = 1000;

%% Plotting the location of the electrodes in 3D
plot3([chanlocs.X], [chanlocs.Y], [chanlocs.Z], 'ko', 'MarkerFaceColor','k');

%Labelling the electrodes
hold on
for tmp = 1:length(electrodes)
    text(chanlocs(tmp).X + 3, chanlocs(tmp).Y, chanlocs(tmp).Z, chanlocs(tmp).labels)
end

xlabel('X'), ylabel('Y'), zlabel('Z')
title('Electrode locations')
axis square

%% Plotting the location of the electrodes in 2D

%TODO: ask if it's ok to use ext sources
%topoplotIndie(, chanlocs, 'electrodes', 'numbers');

title('2D topographical map');

%% PSD

%Plotting the differences in PSD between working and resting conditions for
%each subject
for signal = (length(myFiles) - 1)/2:-1:1 %Same dimension between ds
    
    %Computing the PSD with the Welch method for every electrode in both
    %working and resting conditions
    for electrode = length(electrodes):-1:1

        %Expliciting the signals in time domain
        timeWorkSig = myWorkData(signal).(electrodes(electrode));
        timeRestSig = myRestData(signal).(electrodes(electrode));

        %Trimming the time series to get rid of signal error
        timeWorkSig = timeWorkSig(1:length(timeWorkSig) - sigErr);
        timeRestSig = timeRestSig(1:length(timeRestSig) - sigErr);

        %Computing the PSD using the Welch method (specifics in references)
        freqWorkSig(electrode, :) = pwelch(timeWorkSig, hamming(Fs*10), Fs*0.1, Fs);
        freqRestSig(electrode, :) = pwelch(timeRestSig, hamming(Fs*10), Fs*0.1, Fs);

    end

    %TODO: decide how to sort the electrodes
    [~,sortXidx] = sort([chanlocs.X]);

    figure

    %Plotting the PSDs for working conditions
    subplot(2, 1, 1)
    imagesc(hz,[], freqWorkSig(sortXidx, 1:length(hz)));
    set(gca, 'xlim', [0 70], 'clim', [0 10]); %Focusing on the 0-70Hz range
    colorbar

    xlabel('Frequency (Hz)'), ylabel('Electrode')
    title(strcat('Subject ', int2str(signal), ' work'))

    %Plotting the PSDs for resting conditions
    subplot(2, 1, 2)
    imagesc(hz,[], freqRestSig(sortXidx, 1:length(hz)));
    set(gca, 'xlim', [0 70], 'clim', [0 10]); %Focusing on the 0-70Hz range
    colorbar

    xlabel('Frequency (Hz)'), ylabel('Electrode')
    title(strcat('Subject ', int2str(signal), ' rest'))

    %TODO: choose colorbar range relatively to the normalization of the
    %fourier transform(?)
        
    %TODO: here plotting the actual differences once having calculated all
    %the psds

end

%% Detrended Fluctuation Analysis

for signal = length(myRestData):-1:1 %Same dimension between ds
    figure;
    for electrode = length(electrodes):-1:1

        %Expliciting the signals in time domain
        timeWorkSig = myWorkData(signal).(electrodes(electrode));
        timeRestSig = myRestData(signal).(electrodes(electrode));

        %Summing
        %TODO: detrending the signal and summing
    end
end