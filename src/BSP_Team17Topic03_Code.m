clear variables; close all; clc

%Setting the path to the resources depending on the current OS
if isunix
    path = "/Users/mattiapezzano/Documents/GitHub/proj-bsp-2023/Data/";
else
    path = "Data\";
end

%Loading the data files into a dir struct
myFiles = dir(strcat(path, "*.mat"));

%Loading the file containing the specifics about the EEG
EEG = load(strcat(path, "chanlocs.mat"));

%Iterating through the activity data (skipping the first element)
for i = (length(myFiles) - 1)/2:-1:1
    filePath = strcat(path, myFiles(i * 2 - 1).name); %Odd spacing w/ given data
    myWorkData(i) = load(filePath);
end

%Iterating through the rest data
for j = (length(myFiles) - 1)/2:-1:1
    filePath = strcat(path, myFiles(j * 2).name); %Even spacing w/ given data
    myRestData(j) = load(filePath);
end

%Expliciting the struct attributes in order to iterate through them
%TODO: see if chanlocs order is the same as the signal one, use chanlocs in
%case it is
electrodes = ["C3","C4","CZ","F3","F4","F7","F8","FP1","FP2","FZ","O1","O2","P3","P4","PZ","T3","T4","T5","T6"];

Fs = 500;
% t = (1:30*Fs)/Fs;

%% Plotting the location of the electrodes in 3D
plot3([EEG.chanlocs.X], [EEG.chanlocs.Y], [EEG.chanlocs.Z], 'ko', 'MarkerFaceColor','k');

hold on
for i = 1:length([EEG.chanlocs.X])
    %TODO: correct locations
    text(EEG.chanlocs(i).X + 3, EEG.chanlocs(i).Y, EEG.chanlocs(i).Z, EEG.chanlocs(i).labels)
end

xlabel('X'), ylabel('Y'), zlabel('Z')
title('Electrode locations')
axis square

%% Plotting the location of the electrodes in 2D


% %TODO:  Undefined function 'finputcheck' for input arguments of type 'cell'.
% topoplot(zeros(length(electrodes), 1), EEG.chanlocs, 'electrodes', 'numbers')

% for electrode = length(electrodes):-1:1
%     figure
%     [handle,Zi,grid,Xi,Yi] = topoplot([chanLocs.chanlocs(i).X], [chanLocs.chanlocs(i).Y], chanLocs.chanlocs(i));
%     plot(grid);
% end

%% PSD

% Plotting the differences in the PSD fucntion between the rest and 
% activity state for each person for each electrode
for signal = (length(myFiles) - 1)/2:-1:1 %Same dimension between ds
    figure;
    for electrode = length(electrodes):-1:1

        %Expliciting the signals in time domain
        timeWorkSig = myWorkData(signal).(electrodes(electrode));
        timeRestSig = myRestData(signal).(electrodes(electrode));

        %Expliciting the signals in the frequency domain using the Welch
        %method as specified in the reference paper
        freqWorkSig = pwelch(timeWorkSig, hamming(Fs*10), Fs*0.1, Fs);
        freqRestSig = pwelch(timeRestSig, hamming(Fs*10), Fs*0.1, Fs);

        %Plotting the overlap between the two power spectrum densities
        subplot(5, 4, electrode);
        plot(freqWorkSig); xlim([0 70]);
        hold on
        plot(freqRestSig); xlim([0 70]);

        %TODO: fix title: keep info about name of file
        title(string(electrodes(electrode)));
    end
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