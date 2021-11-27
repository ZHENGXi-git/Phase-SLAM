classdef PoseNode < handle
    %POSENODE A class for pose node
    
    properties (Access = public)
        id    % Id of this pose node
        pose  % Pose of this pose node
    end  % properties public
    
    properties (Dependent = true)
        x    % X coordinate
        y    % Y coordinate
        z    % Z coordinate
        alpha 
        beta  % Yaw angle
        gamma
        Rt   % Transformation local to global
    end  % properties dependent
    
    methods
        
        function obj = PoseNode(id, pose)
            % Constructor of PoseNode
            obj.id   = id;
            obj.pose = pose(:);
        end
        
        function plot(obj)
            % Plot all pose nodes position
            x = [obj.x];
            y = [obj.y];
            z = [obj.z];
            plot(x, y, z);
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

        function beta = get.beta(obj)
            beta = obj.pose(4);
        end
        function gamma = get.gamma(obj)
            gamma = obj.pose(4);
        end
        
        function Rt = get.Rt(obj)    
            cos_alpha = cos(obj.alpha);
            sin_alpha = sin(obj.alpha);
            cos_beta = cos(obj.beta);
            sin_beta = sin(obj.beta);
            cos_gamma = cos(obj.gamma);
            sin_gamma = sin(obj.gamma);
            
            R11 = cos_beta * cos_gamma;
            R12 = -sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma;
            R13 = -cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma;
            R21 = cos_beta * sin_gamma;
            R22 = -sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma;
            R23 = -cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma;
            R31 = sin_beta;
            R32 = sin_alpha * cos_beta;
            R33 = cos_alpha * cos_beta;
   
            Rt =[R11  R12  R13  v(1)   % maybe change!!
                 R21  R22  R23  v(2)
                 R31  R32  R33  v(3)
                  0    0    0    1  ];

        end
        
    end  % methods public
    
end  % classdef
