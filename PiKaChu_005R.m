%% Real Dataset-PiKaChu005

close all; clear;
addpath(genpath('.\lib'));

LoadRoad = 'dataset/Real_dataset/PiKaChu/';
ShowResult = 1;  % show the reconstruction result ?

PhaseSLAM_RobotArm(LoadRoad, ShowResult);