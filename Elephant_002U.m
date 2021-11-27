%% Unreal Dataset-elephant_001

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Unreal_dataset/elephant/';

ShowResult = 1;   % show the reconstruction result ?

PhaseSLAM_Unreal(LoadRoad, ShowResult);