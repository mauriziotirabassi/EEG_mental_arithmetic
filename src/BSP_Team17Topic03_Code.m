clear variables; close all; clc

%Setting the path to the resources relatively to the running OS
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

subjNum = (length(myFiles) - 1)/2; % Number of subjects ("- 1" to exclude chanlocs.mat)

%Iterating through the activity data (skipping the first element)
for tmp = subjNum:-1:1 
    filePath = strcat(path, myFiles(tmp * 2 - 1).name); %Odd spacing w/ given data
    myWorkData(tmp) = load(filePath);
end

%Iterating through the rest data
for tmp = subjNum:-1:1
    filePath = strcat(path, myFiles(tmp * 2).name); %Even spacing w/ given data
    myRestData(tmp) = load(filePath);
end

chanNum = 19; %Number of electrodes
Fs = 500; %Sampling frequency
% hz = linspace(0, Fs/2, Fs*10); %Frequency scale
sigErr = 1000; %Number of samples to ignore
dimW = 45; %Window length in seconds to approximate for stationarity

freqBands = [1, 4, 8, 13, 20, 30, 40]; %Relevant frequency bands

%Defining the vector of the electrode channels
electrodes = strings(1, chanNum);
for i = 1:chanNum %TODO: length chanlocs - 2 redundant
    electrodes(i) = upper(chanlocs(i).labels);
end

workBands = zeros(length(electrodes), length(freqBands), subjNum);
restBands = zeros(length(electrodes), length(freqBands), subjNum);

%% Plotting the location of the electrodes in 3D
figure
plot3([chanlocs.X], [chanlocs.Y], [chanlocs.Z], 'ko', 'MarkerFaceColor','k');

%Labelling the electrodes
hold on
for tmp = 1:length(electrodes)
    text(chanlocs(tmp).X + 3, chanlocs(tmp).Y, chanlocs(tmp).Z, chanlocs(tmp).labels)
end

xlabel('X'), ylabel('Y'), zlabel('Z')
title('Electrode Positions')
axis square

%% PSD

%Plotting the differences between the power spectral densities in work and
%rest conditions for each subject
figure
for subj = subjNum:-1:1
    
    %Computing the PSD with the Welch method for every electrode in both
    %working and resting conditions
    for electrode = 1:length(electrodes)

        %Expliciting the signals in time domain
        timeWorkSig = myWorkData(subj).(electrodes(electrode));
        timeRestSig = myRestData(subj).(electrodes(electrode));

        %Trimming the time series to get rid of signal error
        timeWorkSig = timeWorkSig(1:length(timeWorkSig) - sigErr);
        timeRestSig = timeRestSig(1:length(timeRestSig) - sigErr);

        %Selecting a limited window to approximate for stationarity
        % (dimW seconds centered in the signal median)
        timeWorkSig = timeWorkSig(length(timeWorkSig)/2 - Fs*(dimW/2) : length(timeWorkSig)/2 + Fs*(dimW/2));
        timeRestSig = timeRestSig(length(timeRestSig)/2 - Fs*(dimW/2) : length(timeRestSig)/2 + Fs*(dimW/2));

        %Computing the PSD using the Welch method (specifics in references)
        [freqWorkSig, hz] = pwelch(timeWorkSig, hamming(Fs*10), Fs*0.1, Fs*10);
        [freqRestSig, ~] = pwelch(timeRestSig, hamming(Fs*10), Fs*0.1, Fs*10);
        
        %Normalizing the psds to highlight the differences between work
        %and rest conditions
        ratioSig(electrode, :, subj) = 10*log10(freqWorkSig./freqRestSig);

        %Averaging over bands of interest:
        % delta(1-4)Hz, theta(4-8)Hz, alpha(8-13)Hz, 
        % beta1(13-20)Hz, beta2(20-30)Hz, gamma(30-40)Hz
        for freq = 1:(length(freqBands) - 1)
            startBand = freqBands(freq);
            endBand = freqBands(freq + 1);
            interval = hz(startBand*10 + 1) : 0.1 : hz(endBand*10 + 1);
            intervalX = linspace(startBand, endBand, (endBand - startBand)*10 + 1);
            intervalY = ratioSig(electrode, startBand*10 + 1: endBand*10 + 1, subj);
            bandSig(electrode, freq, subj) = trapz(intervalX, intervalY);
        end

        %Definig the double ds
        workRecord(electrode, :, subj) = timeWorkSig;
        restRecord(electrode, :, subj) = timeRestSig;

    end

    %Denormalizing the frequency range (implicitly normalized by the pwelch
    %function)
    hz = hz*Fs/2/pi;

    %Sorting the electrodes relatively to their X coordinate
    % (from the occipital lobe to the frontal lobe)
    [~, sortXidx] = sort([chanlocs.X]);

    %Plotting the ratio normalized psd
    subplot(3, 2, subj)
    imagesc(hz, [], ratioSig(sortXidx, 1:length(hz), subj));
    set(gca, 'xlim', [0 70], 'clim', [-30 30]); %Focusing on the 0-70Hz range
    colorbar
    xlabel('Frequency (Hz)'), ylabel('F <-- --> O')
    title(strcat('Subj', int2str(subj), ' Work'))
    
    %TODO: choose colorbar range relatively to the normalization of the
    %fourier transform(?): run a simulation to know the actual value
    %beforehand
end

%% Plotting the average topographical maps
figure
avgBandSig = mean(bandSig, 3);
for freq = 1:(length(freqBands) - 1)
    subplot(2, 3, freq)
    topoplot(avgBandSig(:, freq), chanlocs)
    colormap parula
    title(['(' int2str(freqBands(freq)) '-' int2str(freqBands(freq + 1)) ')Hz'])
end

%% Averaging operations

%TODO: average of the PSD for the 4 given ranges. Four topoplots for each
%patient both in rest and working conditions

%Dividing the PSDs relatively to the lobe they belong to
for subj = subjNum:-1:1
    frontal = 7;
    temporalDx = 2;
    temporalSx = 2;
    parietal = 6;
    occipital = 2;

    for electrode = 1:length(electrodes)
        eTag = electrodes(electrode);

        if(startsWith(eTag, 'F'))
            frontalWork(frontal, :, subj) = freqWorkSig(electrode, :, subj);
            frontalRest(frontal, :, subj) = freqRestSig(electrode, :, subj);
            frontal = frontal - 1;
        elseif(startsWith(eTag, 'C') || startsWith(eTag, 'P'))
            parietalWork(parietal, :, subj) = freqWorkSig(electrode, :, subj);
            parietalRest(parietal, :, subj) = freqRestSig(electrode, :, subj);
            parietal = parietal - 1;
        elseif(startsWith(eTag, 'O'))
            occipitalWork(occipital, :, subj) = freqWorkSig(electrode, :, subj);
            occipitalRest(occipital, :, subj) = freqRestSig(electrode, :, subj);
            occipital = occipital - 1;
        elseif(eTag == "T4" || eTag == "T6")
            temporalDxWork(temporalDx, :, subj) = freqWorkSig(electrode, :, subj);
            temporalDxRest(temporalDx, :, subj) = freqRestSig(electrode, :, subj);
            temporalDx = temporalDx - 1;
        elseif(eTag == "T3" || eTag == "T5")
            temporalSxWork(temporalSx, :, subj) = freqWorkSig(electrode, :, subj);
            temporalSxRest(temporalSx, :, subj) = freqRestSig(electrode, :, subj);
            temporalSx = temporalSx - 1;
        end
    end

    %Computing the lobe mean
    fWMean(subj, :) = mean(frontalWork(:, :, subj));
    fRMean(subj, :) = mean(frontalRest(:, :, subj));

    pWMean(subj, :) = mean(parietalWork(:, :, subj));
    pRMean(subj, :) = mean(parietalRest(:, :, subj));

    oWMean(subj, :) = mean(occipitalWork(:, :, subj));
    oRMean(subj, :) = mean(occipitalRest(:, :, subj));

    tDWMean(subj, :) = mean(temporalDxWork(:, :, subj));
    tDRMean(subj, :) = mean(temporalDxRest(:, :, subj));

    tSWMean(subj, :) = mean(temporalSxWork(:, :, subj));
    tSRMean(subj, :) = mean(temporalSxRest(:, :, subj));

    %Plotting
    figure
    plot(hz, fWMean(subj, :))
    hold on
    plot(hz, fRMean(subj, :))

end

%% 
%Dividere freqWorkSig per freqRestSig (normalizzazione rispetto a baseline)

%Differenze tra rest e sig tra media di tutti subject (implicita se
%normalizzato)

%

%% topoplots

%% coherence between channels
for k=subjNum:-1:1
    for i= 1:length(electrodes)
        for j= 1:length(electrodes)
            RestMSC(i,j,:,k)=mscohere(RestRecord(i,:,k),RestRecord(j,:,k),hamming(Fs*10), Fs*0.1, Fs*10);
            WorkMSC(i,j,:,k)=mscohere(WorkRecord(i,:,k),WorkRecord(j,:,k),hamming(Fs*10), Fs*0.1, Fs*10);
        end
    end
end

%Filtro e poi psd o dividi a occhio la finestra?
