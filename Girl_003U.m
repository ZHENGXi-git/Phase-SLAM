%% Unreal Dataset-girl_001

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Unreal_dataset/girl/';
ShowResult = 1;   % show the reconstruction result ?

PhaseSLAM_Unreal(LoadRoad, ShowResult);