classdef PoseNode_v < handle
    %POSENODE A class for pose node
    
    properties (Access = public)
        id    % Id of this pose node
        pose  % Pose of this pose node
    end  % properties public
    
    properties (Dependent = true)
        x    % X coordinate
        y    % Y coordinate
        z    % Z coordinate
        alpha  % Alpha angle
        rt   % Transformation local to global
    end  % properties dependent
    
    methods
        
        function obj = PoseNode_v(id, pose)
            % Constructor of PoseNode
            obj.id   = id;
            obj.pose = pose(:);
        end
        
        function plot(obj)
            % Plot all pose nodes position
            x = [obj.x];
            y = [obj.y];
            z = [obj.z];
            plot(y, z);
        end
        
        function x = get.x(obj)
            x = obj.pose(1);
        end
        
        function y = get.y(obj)
            y = obj.pose(2);
        end
        
        function z = get.z(obj)
            z = obj.pose(3);
        end
        
        function alpha = get.alpha(obj)
            alpha = obj.pose(4);
        end
        
        function rt = get.rt(obj)    
            R = [1        0             0
                 0   cos(obj.alpha)  sin(obj.alpha)
                 0  -sin(obj.alpha)  cos(obj.alpha)];
            rt = [R [obj.x; obj.y; obj.z];0 0 0 1];
        end
        
    end  % methods public
    
end  % classdef
