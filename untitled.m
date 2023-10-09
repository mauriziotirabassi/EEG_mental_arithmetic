clear all

%TODO: implement matrices ds

%Importing the rest data
for i = 1:6
    s(i) = load(strcat('Documentation\Data\Subject0', num2str(i), '_1.mat'));
end

%Importing the activity data
for i = 1:6
    t(i) = load(strcat('Documentation\Data\Subject0', num2str(i), '_2.mat'));
end

