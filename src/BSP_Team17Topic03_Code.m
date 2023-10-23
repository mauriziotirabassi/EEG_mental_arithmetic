clear variables; close all; clc

% Setting the path to the resources relatively to the running OS
if isunix
    path = '/Users/mattiapezzano/Documents/GitHub/proj-bsp-2023/src/Data/';
else
    path = 'Data\';
end

myFiles = dir(strcat(path, '*.mat')); % Loading the data files into a dir struct
load(strcat(path, 'chanlocs.mat')); % Loading the file containing the specifics about the EEG

subjNum = (length(myFiles) - 1)/2; % Number of subjects ("- 1" to exclude chanlocs.mat)

% Iterating through the activity data (skipping the first element)
for tmp = subjNum:-1:1 
    filePath = strcat(path, myFiles(tmp * 2 - 1).name); % Odd spacing w/ given data
    myWorkData(tmp) = load(filePath);
end

% Iterating through the rest data
for tmp = subjNum:-1:1
    filePath = strcat(path, myFiles(tmp * 2).name); % Even spacing w/ given data
    myRestData(tmp) = load(filePath);
end

chanNum = 19; % Number of electrodes
Fs = 500; % Sampling frequency
sigErr = 1000; % Number of samples to ignore
dimW = 45; % Window length in seconds to approximate for stationarity

freqBands = [1, 4, 8, 13, 20, 30, 40]; % Relevant frequency bands

% Defining the vector of the electrode channels
electrodes = strings(1, chanNum);
for i = 1:chanNum %TODO: length chanlocs - 2 redundant
    electrodes(i) = upper(chanlocs(i).labels);
end

%% PLOTTING 3D ELECTRODE LOCATION
figure
plot3([chanlocs.X], [chanlocs.Y], [chanlocs.Z], 'ko', 'MarkerFaceColor','k');

% Labelling the electrodes
hold on
for tmp = 1:length(electrodes)
    text(chanlocs(tmp).X + 3, chanlocs(tmp).Y, chanlocs(tmp).Z, chanlocs(tmp).labels)
end

xlabel('X'), ylabel('Y'), zlabel('Z')
title('Electrode Positions')
axis square

%% DATA PREPROCESSING

% Preallocating memory for the data structures containing the time series
% for each electrode for each subject
workBands = zeros(length(electrodes), length(freqBands), subjNum);
restBands = zeros(length(electrodes), length(freqBands), subjNum);

figure
for subj = subjNum:-1:1

    % Defining the number of electrode channels belonging to each lobe
    frntNum = 7; prtNum = 6; occNum = 2; tmpdxNum = 2; tmpsxNum = 2;

    for electrode = 1:length(electrodes)
        eTag = electrodes(electrode); % Electrode channel tag

        % Expliciting the signals in time domain
        timeWork = myWorkData(subj).(eTag);
        timeRest = myRestData(subj).(eTag);

        % Trimming the time series to get rid of signal error
        timeWork = timeWork(1:length(timeWork) - sigErr);
        timeRest = timeRest(1:length(timeRest) - sigErr);

        % Selecting a limited window to approximate for stationarity
        % (dimW seconds centered in the signal median)
        timeWork = timeWork(length(timeWork)/2 - Fs*(dimW/2) : length(timeWork)/2 + Fs*(dimW/2));
        timeRest = timeRest(length(timeRest)/2 - Fs*(dimW/2) : length(timeRest)/2 + Fs*(dimW/2));

        % Saving the time signals into a data structure
        workRecord(electrode, :, subj) = timeWork;
        restRecord(electrode, :, subj) = timeRest;

        % Computing the PSD using the Welch method (specifics in references)
        [psdWork, hz] = pwelch(timeWork, hamming(Fs*10), Fs*0.1, Fs*10);
        [psdRest, ~] = pwelch(timeRest, hamming(Fs*10), Fs*0.1, Fs*10);
        
        % Normalizing the psds to highlight the differences between work
        % and rest conditions
        psdRatio(electrode, :, subj) = 10*log10(psdWork./psdRest);

        % Averaging over bands of interest:
        % delta(1-4)Hz, theta(4-8)Hz, alpha(8-13)Hz, 
        % beta1(13-20)Hz, beta2(20-30)Hz, gamma(30-40)Hz
        for freq = 1:(length(freqBands) - 1)

            % Defining the frequency bands
            startBand = freqBands(freq);
            endBand = freqBands(freq + 1);

            % Defining the integration intervals
            intervalX = linspace(startBand, endBand, (endBand - startBand)*10 + 1);
            intervalY = psdRatio(electrode, startBand*10 + 1: endBand*10 + 1, subj);

            % Computing the mean frequency for every band for every subject
            psdBand(electrode, freq, subj) = trapz(intervalX, intervalY);

        end

        %TODO: decide if it's more convenient having it isolated inside
        %another section

        % Dividing electrode channels among lobes for each subject
        if(startsWith(eTag, 'F'))
            frntPsd(frntNum, :, subj) = psdRatio(electrode, :, subj);
            frntNum = frntNum - 1;
        elseif(startsWith(eTag, 'P') || startsWith(eTag, 'C'))
            prtPsd(prtNum, :, subj) = psdRatio(electrode, :, subj);
            prtNum = prtNum - 1;
        elseif(startsWith(eTag, 'O'))
            occPsd(occNum, :, subj) = psdRatio(electrode, :, subj);
            occNum = occNum - 1;
        elseif(eTag == "T4" || eTag == "T6")
            tmpdxPsd(tmpdxNum, :, subj) = psdRatio(electrode, :, subj);
            tmpdxNum = tmpdxNum - 1;
        elseif(eTag == "T3" || eTag == "T5")
            tmpsxPsd(tmpsxNum, :, subj) = psdRatio(electrode, :, subj);
            tmpsxNum = tmpsxNum - 1;
        end

    end

    % Computing the average frequency for each lobe for each subject
    frntPsdAvg(subj, :) = mean(frntPsd(:, :, subj));
    prtPsdAvg(subj, :) = mean(prtPsd(:, :, subj));
    occPsdAvg(subj, :) = mean(occPsd(:, :, subj));
    tmpdxPsdAvg(subj, :) = mean(tmpdxPsd(:, :, subj));
    tmpsxPsdAvg(subj, :) = mean(tmpsxPsd(:, :, subj));

    % Denormalizing the frequency range (implicitly normalized by the pwelch
    % function)
    %TODO: if not plotting can define outside of the outer for loop
    hz = hz*Fs/2/pi;

    %TODO: decide if graphs are relevant (include process decision making
    %in the paper)

    % Sorting the electrodes relatively to their X coordinate
    % (from the occipital lobe to the frontal lobe)
    [~, sortXidx] = sort([chanlocs.X]);

    % Plotting the ratio normalized psd
    subplot(3, 2, subj)
    imagesc(hz, [], psdRatio(sortXidx, 1:length(hz), subj));

    % Deciding the color range knowing a priori the frequency range
    set(gca, 'xlim', [0 70], 'clim', [-30 30]); %Focusing on the 0-70Hz range
    colorbar
    xlabel('Frequency (Hz)'), ylabel('F <-- --> O')
    title(strcat('Subj', int2str(subj), ' Work'))
    
end

%% PLOTTING AVERAGED TOPOGRAPHICAL MAPS

%Computing the average over subjects
avgBandSig = mean(psdBand, 3);

figure
for freq = 1:(length(freqBands) - 1)
    subplot(2, 3, freq)
    topoplot(avgBandSig(:, freq), chanlocs)
    colormap parula
    title(['(' int2str(freqBands(freq)) '-' int2str(freqBands(freq + 1)) ')Hz'])
end

%% COHERENCE BETWEEN CHANNELS

for k=subjNum:-1:1
    for i= 1:length(electrodes)
        for j= 1:length(electrodes)
            RestMSC(i,j,:,k)=mscohere(RestRecord(i,:,k),RestRecord(j,:,k),hamming(Fs*10), Fs*0.1, Fs*10);
            WorkMSC(i,j,:,k)=mscohere(WorkRecord(i,:,k),WorkRecord(j,:,k),hamming(Fs*10), Fs*0.1, Fs*10);
        end
    end
end
