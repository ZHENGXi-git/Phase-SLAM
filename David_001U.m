%% Unreal Dataset-david001

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Unreal_dataset/david/';

ShowResult = 1;  % show the reconstruction result ?

PhaseSLAM_Unreal(LoadRoad, ShowResult);
