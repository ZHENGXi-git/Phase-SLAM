%% Real Dataset-david004

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Real_dataset/David/';
ShowResult = 1;  % show the reconstruction result ?

PhaseSLAM_RobotArm(LoadRoad, ShowResult);