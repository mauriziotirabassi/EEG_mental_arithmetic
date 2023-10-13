clear variables; close all; clc

Fs = 500;
% t = (1:30*Fs)/Fs;

%TODO: when using MACOS explicit the full length path: 
%path="/Users/mattiapezzano/Documents/GitHub/proj-bsp-2023/Data/";
path="Data\";
myFiles = dir(strcat(path, "*.mat"));

%TODO: fix path from directory
chanLocs = load(strcat(path, "chanlocs.mat"));

%Iterating through the activity data (skipping the first element)
for i = 6:-1:1
    filePath = strcat(path, myFiles(i * 2 - 1).name); %Odd spacing w/ given data
    myWorkData(i) = load(filePath);
end

%Iterating through the rest data
for j = 6:-1:1
    filePath = strcat(path, myFiles(j * 2).name); %Even spacing w/ given data
    myRestData(j) = load(filePath);
end

%Expliciting the struct attributes in order to iterate through them
electrodes = ["C3","C4","CZ","F3","F4","F7","F8","FP1","FP2","FZ","O1","O2","P3","P4","PZ","T3","T4","T5","T6"];

%% Plotting the location of the electrodes in 3D
plot3([chanLocs.chanlocs.X], [chanLocs.chanlocs.Y], [chanLocs.chanlocs.Z], 'ko', 'MarkerFaceColor','k');

hold on
for i = 1:length([chanLocs.chanlocs.X])
    %TODO: correct locations
    text(chanLocs.chanlocs(i).X + 3, chanLocs.chanlocs(i).Y, chanLocs.chanlocs(i).Z, chanLocs.chanlocs(i).labels)
end

xlabel('X'), ylabel('Y'), zlabel('Z')
title('Electrode locations')
axis square

%% Plotting the location of the electrodes in 2D
%TODO: topoplot

%% PSD

% Plotting the differences in the PSD fucntion between the rest and 
% activity state for each person for each electrode
for signal = length(myRestData):-1:1 %Same dimension between ds
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

        %TODO: fix title
        title(strcat(int2str(signal), ":", electrodes(electrode)));
    end
end

