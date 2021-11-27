%% Real Dataset-Sona006

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Real_dataset/Sona/';
ShowResult = 1;  % show the reconstruction result ?

PhaseSLAM_RobotArm(LoadRoad, ShowResult);