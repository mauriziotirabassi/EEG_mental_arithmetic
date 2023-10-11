clear variables; close all; clc

% Fs = 500;
% t = (1:30*Fs)/Fs;

%Iteration from the last element is relevant for dynamic memory allocation

%TODO: when using MACOS explicit the full length path: 
% /Users/mattiapezzano/Documents/GitHub/proj-bsp-2023/Documentation/Data/
myFiles = dir("Documentation\Data\*.mat");

%Ierating through the activity data
for i = (length(myFiles) - 1):-2:1 %Odd spacing w/ given data
    filePath = strcat("Documentation\Data\", myFiles(i).name);
    myActivityData(i) = load(filePath);
end

%Ierating through the rest data
for i = length(myFiles):-2:1 %Even spacing w/ given data
    filePath = strcat("Documentation\Data\", myFiles(i).name);
    myRestData(i) = load(filePath);
end

%Expliciting the struct attributes in order to iterate through them
electrodes = ["C3","C4","CZ","F3","F4","F7","F8","FP1","FP2","FZ","O1","O2","P3","P4","PZ","T3","T4","T5","T6"];

%Plotting the differences(?) between the rest and activity state for each
%person for each electrode
for signal = length(myActivityData):-1:1 %Same dimension between ds
    for electrode = length(electrodes):-1:1
        data = myActivityData(signal).(electrodes(electrode));
        %TODO:superimposition of plots not feasible due to different length
        figure;
        plot(data);
        title(strcat(int2str(signal), ":", electrodes(electrode)));
    end
end


%%SAMPLE STRUCT PARAM IT
% fn=fieldnames(structure);
% %loop through the fields
% for i=1: numel(fn)
%     fn1=fieldnames(structure.(fn{i}))
%     for j=1: numel(fn1)
%         fn2=fieldnames(structure.(fn{i}).(fn1{j}));
%         for k=1: numel(fn2)
%             %access the data
%             data=structure.(fn{i}).(fn1{j}).(fn2{k});     
%         end 
%     end
% end

%%SAMPLE SIG PROC
% clear;clc;
% % Generate white gaussian noise
% time_req = 10; % seconds
% fs = 100; %Hz
% samples = time_req*fs; % seconds
% new_time = linspace(0,time_req,samples);
% noise = randn(1,samples); % AWGN
% b=1;
% a=[1 -0.3 0.7];
% filt_noise = filter(b,a,noise);
% 
% % Noise and PSD
% figure; subplot(2,2,1)
% plot(new_time,noise); xlabel('samples'); ylabel('amplitude');
% title('White Gaussian Noise')
% subplot(2,2,2)
% periodogram(noise); xlabel('Normalized Frequency'); ylabel('Power');
% title('PSD of White Gaussian Noise')
% subplot(2,2,3)
% plot(new_time,filt_noise); xlabel('samples'); ylabel('amplitude');
% title('Output of the system')
% subplot(2,2,4)
% periodogram(filt_noise);  xlabel('Normalized Frequency'); ylabel('Power');
% title('PSD of the output')
