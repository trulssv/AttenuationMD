clear all, close all, clc

%% Path to data

% Endoleak contrast no leak

dataroot = '/mnt/data/MDdata/Images_Endoleak/20260224-e559/s310/';
file = 'i000300.dcm';
datapath = strcat(dataroot, file);

%% Load data

i = dicominfo(datapath);

% For whatever reason the MDs are stored here...

i1 = double(reshape(i.Private_0053_1095, [512, 512]));
i2 = double(reshape(i.Private_0053_1096, [512, 512]));
i3 = double(reshape(i.Private_0053_1097, [512, 512]));
%% Some data processing

i1(i1>5e3) = 0;
i2(i2>5e3) = 0;
i3(i3>5e3) = 0;

% 1
% 0.1 mg/ml for 0.2 mg/ml for I