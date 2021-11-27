function[] = PhaseSLAM_RobotArm(LoadRoad, SHOWRESULT)

% camera matrix
row = 1024;
col = 1280;
fc1 = 3355.302388053488812;
fc2 = 3355.064851430768158;
cc1 = 650.807692865668059;
cc2 = 524.954727398199907;
alphac = -0.000028723621966;
fp1 = 2193.252320427087398;
fp2 = 2192.481018734629288;
cp1 = 620.776653065079813;
cp2 = 523.502024084542313;
alphap = -0.001396388086297;
Inc = [fc1  alphac*fc1  cc1;
       0       fc2      cc2;
       0        0       1];
Rc = [ 0.986554368359845, -0.001449007437864633,   -0.1634269826177788
      0.000943781303290944,     0.999994533190938, -0.003169047978393839
      0.1634306811677622,  0.002972198795978773,     0.986550342601612];
Tc = [80.13418712805978
     -0.9181180126974828
      37.86801174501727];
% projector matrix
Rp = eye(3); 
Tp = [0 0 0]';
Inp = [fp1  alphap*fp1  cp1
       0       fp2      cp2
       0        0       1];

CameraM = Inc * [Rc, Tc];
ProjectorM = Inp * [Rp, Tp];

%% Initialization
InputData = load([LoadRoad, 'InputData.txt']);
PoseNodeNum = length(InputData);
nodes = zeros(4, PoseNodeNum + 1);
edges = zeros(5, PoseNodeNum + 1);
Rn = eye(3);
tn = [0 0 0]';
R_t.R = Rn;
R_t.t = tn;
count = 1;
Interval = 1;
Length = row * col;
GaussM = randn([100 Length]);
threshold = 100;
KeyFrame = 1;

%% Phase Image 
Phase1 = load([LoadRoad, 'PhaseData1.txt']);

PhaseLine = reshape(Phase1, [], 1);
CompressPhase = GaussM * PhaseLine;

[X1, Y1, Z1] = Compute3D(Phase1, CameraM, ProjectorM, col, row, 0);

TotalTime = 0.0;
LoopTime = 0.0;
figure('name','Reconstruction without Loop Closure');
for i = 1: 1: PoseNodeNum

    if(SHOWRESULT == 1)
        [X, Y, Z] = Compute3D(Phase1, CameraM, ProjectorM, col, row, 2);
        PointCloud.X = X;
        PointCloud.Y = Y;
        PointCloud.Z = Z;
        PointList(i) = PointCloud;
    end
    idx = InputData(i, 1);   
    Phase2 = load([LoadRoad, 'PhaseData', num2str(idx + 1),'.txt']);  
    [X2, Y2, Z2] = Compute3D(Phase2, CameraM, ProjectorM, col, row, 0);

    ROI = InputData(i, (2: 5));
    initialization = InputData(i, (6: 17));
    iteration = 5000;
    t = clock;
    [motion, diff1] = OneOrderOptimizer(Phase1, Phase2, X1, Y1, Z1, ROI, Inp, CameraM, initialization, iteration, 1);
    TotalTime = etime(clock, t) + TotalTime;
    fprintf('The %dth pose estimation cost: %d\n', i, diff1);
    [pose, Rn, tn, R_t] = PoseCompution(motion, Rn, tn, R_t, count, 1);
    count = count + 1;
    nodes(:, i + 1) = [i, pose];
    edges(:, i) = [i, i - 1, motion(2), motion(3), motion(4)];
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
                edges(:, i + 1) = [0, i, 0, 0, 0];
                pg = PoseGraph_v();
                pg.readGraph_v(nodes, edges);
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
        [X2, Y2, Z2] = Compute3D(Phase2, CameraM, ProjectorM, col, row, 2);
        [X2, Y2, Z2] = PCFusion(X2, Y2, Z2, R_t, count-1, 2);
        plot3(X2(:),Y2(:),Z2(:),'.b','MarkerSize',1); grid on; axis equal;
        xlim([-150 350]); ylim([-200 200]);zlim([300 800]);
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
        newR_t = v2m_v(pg.pose(:, i));
        [NewX, NewY, NewZ] = PCFusion(X, Y, Z,newR_t, 1, 3);
        PointCloud.X = NewX;
        PointCloud.Y = NewY;
        PointCloud.Z = NewZ;

        plot3(NewX(:), NewY(:), NewZ(:),'.b','MarkerSize',1); grid on; axis equal;
        xlim([-150 350]); ylim([-200 200]);zlim([300 800]);
        title('The Reconstruction Result with Loop')
        hold on;
    end
end

end