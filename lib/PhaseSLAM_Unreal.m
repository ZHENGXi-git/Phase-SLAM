
function[] = PhaseSLAM_Unreal(LoadRoad, SHOWRESULT)

%% The calibration parameters of camera and projector
row = 480;  
col = 640;

Inc = [355.396004745342, 0, 320;    % camera internal parameter
       0, 355.396004745342, 240;    
       0, 0, 1];
Inp = [476.701437037684, 0, 320;    % projector internal parameter
       0, 476.701437037684, 240;
       0, 0, 1];
CameraM = [355.396004745342, 0, 320, 0;
           0, 355.396004745342, 240, 3553.96004745342;
           0, 0, 1, 0];
ProjectorM = [476.701437037684, 0, 320, 0;
              0, 476.701437037684, 240, 0;
              0, 0, 1, 0];

%% Initialization
InputData = load([LoadRoad, 'InputData.txt']);
PoseNodeNum = length(InputData);
nodes = zeros(7, PoseNodeNum + 1);
edges = zeros(8, PoseNodeNum + 1);
Rn = eye(3);
tn = [0 0 0]';
R_t.R = Rn;
R_t.t = tn;
count = 1;
Interval = 1;
Length = row * col;
GaussM = randn([100 Length]);
threshold = 1;
KeyFrame = 1;


%% Phase Image 
Phase1 = load([LoadRoad, 'PhaseData1.txt']);
% Compress sensing
PhaseLine = reshape(Phase1, [], 1);
CompressPhase = GaussM * PhaseLine;

[X1, Y1, Z1] = Compute3D(Phase1, CameraM, ProjectorM, col, row, 1);

TotalTime = 0.0;
LoopTime = 0.0;

figure('name','Reconstruction without Loop Closure');
for i = 1: 1: PoseNodeNum
    
    if(SHOWRESULT == 1)
        [X, Y, Z] = Compute3D(Phase1, CameraM, ProjectorM, col, row, 3);
        PointCloud.X = X;
        PointCloud.Y = Y;
        PointCloud.Z = Z;
        PointList(i) = PointCloud;
    end  
   
    %% compute the phase image and 3D points   
    idx = InputData(i, 1);    
    Phase2 = load([LoadRoad, 'PhaseData', num2str(idx + 1),'.txt']);    
    [X2, Y2, Z2] = Compute3D(Phase2, CameraM, ProjectorM, col, row, 1);  
    
    ROI = InputData(i, (2: 5));
    initialization = InputData(i, (6: 17));
    iteration = 5000;
    t = clock;
    [motion, diff1] = OneOrderOptimizer(Phase1, Phase2, X1, Y1, Z1, ROI, Inp, CameraM, initialization, iteration, 0);
    TotalTime = etime(clock, t) + TotalTime;

    fprintf('The %dth pose estimation cost: %d\n', i, diff1);
    [pose, Rn, tn, R_t] = PoseCompution(motion, Rn, tn, R_t, count, 0);
    count = count + 1;
    nodes(:, i + 1) = [i, pose];
    edges(:, i) = [i, i - 1, motion];
    
    % loop detection
    t0 = clock;
    if (rem(count, Interval) == 0)
        KeyFrame = KeyFrame + 1;
        PhaseLine = reshape(Phase2, [], 1);
        CompressPhase2 = GaussM * PhaseLine;
        CompressPhase(:, KeyFrame) = CompressPhase2;
        for j = 1: 1: KeyFrame-1
            diff = sum(abs(CompressPhase(:, j) - CompressPhase2));
            if (diff < threshold)
                fprintf('The loop closure is detected\n');
                edges(:, i + 1) = [0, i, 0, 0, 0, 0, 0, 0];
                pg = PoseGraph();
                pg.readGraph(nodes, edges);
                pg.lambda = 0;
                pg.optimize(3, true);
                break;
            end
        end
    end
    LoopTime = etime(clock, t0) + LoopTime;
    Phase1 = Phase2;
    X1 = X2;
    Y1 = Y2;
    Z1 = Z2;

    if (SHOWRESULT == 1)
        [X2, Y2, Z2] = Compute3D(Phase2, CameraM, ProjectorM, col, row, 3);
        [X2, Y2, Z2] = PCFusion(X2, Y2, Z2, R_t, count-1, 2);
        plot3(X2(:),Y2(:),Z2(:),'.b','MarkerSize',1); grid on; axis equal;
        xlim([-80 80]); ylim([-60 60]);zlim([80 180]);
        title('The Reconstruction Result without Loop')
        hold on;
    end

end

fprintf('the loop closure time = %d s\n', LoopTime);
fprintf('the Total time = %d s\n', TotalTime);

if (SHOWRESULT == 1)   
    figure('name','Full Reconstruction Result'); 
    
    for i = 1: 1: PoseNodeNum
        X = PointList(i).X;
        Y = PointList(i).Y;
        Z = PointList(i).Z;
        newR_t = v2m(pg.pose(:, i));
        [NewX, NewY, NewZ] = PCFusion(X, Y, Z,newR_t, 1, 1); 
        PointCloud.X = NewX;
        PointCloud.Y = NewY;
        PointCloud.Z = NewZ;

        plot3(NewX(:), NewY(:), NewZ(:),'.b','MarkerSize',1); grid on; axis equal;
        xlim([-80 80]); ylim([-60 60]);zlim([80 180]);
        title('The Reconstruction Result with Loop')
        hold on;

    end
end

end