classdef PoseGraph < handle
    %POSEGRAPH A class for doing pose graph optimization
    
    properties (SetAccess = private)
        node  % Pose nodes in graph
        edge  % Edge in graph
        H     % Information matrix
        b     % Information vector
        E     % Cost Function    
        f     % LM paramter
        q     % LM paramter
        J     %Jacobian
    end  % properties set private
    
    properties (Access = public)
        lambda % LM paramter
    end  % properties public
    
    properties (Dependent = true)
        n_node  % Number of nodes in graph
        n_edge  % Number of edges in graph
        pose    % Poses of all nodes
    end  % properties dependent
    
    methods
        
        function obj = PoseGraph()
            % Constructor of PoseGraph
            obj.node = PoseNode.empty;
            obj.edge = PoseEdge.empty;
        end
        
        function readGraph(obj, nodes, edges)
                  
          %  vertices = vertice.nodes;
            for i_node = 1:size(nodes,2)
                vi = nodes(:,i_node);
                id = vi(1) + 1;
                pose_(1) = vi(2);
                pose_(2) = vi(3);
                pose_(3) = vi(4);
                pose_(4) = vi(5);  
                pose_(5) = vi(6); 
                pose_(6) = vi(7);
                obj.node(i_node) = PoseNode(id, pose_);  % x y yaw
            end

       %     edges = edge.edge;
            for i_edge = 1:size(edges,2)
                ei = edges(:,i_edge);
                id_from   = ei(1) + 1;
                id_to     = ei(2) + 1;
                mean(1)   = ei(3);
                mean(2)   = ei(4);
                mean(3)   = ei(5);
                mean(4)   = ei(6);  % rad
                mean(5)   = ei(7);  % rad
                mean(6)   = ei(8);  % rad
                infm = eye(6,6);
            %    infm
                obj.edge(i_edge) = PoseEdge(id_from, id_to, mean, infm);
            end
            
      %      fprintf('Edges loaded from: %s\n', vfile);
      %      obj.plot();
        end
        
        function plot(obj)
            % Plots pose graph
            hold on;
            obj.node.plot();   % PoseNode.plot();
        end
        
        function optimize(obj, n_iter, vis)
            % Pose graph optimization
            if nargin < 3, vis = false; end
            if nargin < 2, n_iter = 1; end
            
            for i_iter = 1:n_iter
         %       fprintf('Pose Graph Optimization, Iteration %d.\n', i_iter);
                obj.iterate();
         %       fprintf('Iteration %d done.\n', i_iter);
                if vis     % update figure
                %    clf;
              %      obj.plot();
               %     title(sprintf('Iteration %d', i_iter));
                    drawnow;
                end
            end
        end
        
        function iterate(obj)
            % One iteration of pose graph optimization
      %      fprintf('Allocating Workspace.\n');
            % Create new H and b matrices each time
            obj.H = zeros(obj.n_node*6);   % 6n x 6n square matrix
            obj.b = zeros(obj.n_node*6,1); % 6n x 1  column vector
            
      %      fprintf('Linearizing.\n');
            obj.linearize_();
            
      %      fprintf('Solving.\n');
            obj.solve();
        end

        function linearize_(obj)
            
         %   e_plot = zeros(4, obj.n_edge);
        %    E_ij = zeros(1, obj.n_edge);
            f_ = zeros(obj.n_node*6, 1); 
            cost = 0;
            obj.E = 0;
            Jacobian = zeros(obj.n_node*6, obj.n_node*6);
            for i_edge = 1:obj.n_edge
            
                ei = obj.edge(i_edge);   %edge
                % Get edge information
                i_node = ei.id_from;
                j_node = ei.id_to;
                T_z = v2m(ei.mean);
                omega = ei.infm;  
                
                % Get node information
                v_i = obj.node(i_node).pose;
                v_j = obj.node(j_node).pose;
                
                i_ind = id2ind(i_node);    % H_spares
                j_ind = id2ind(j_node);
                
                T_i = v2m(v_i);
                T_j = v2m(v_j);
                R_i = T_i(1:3,1:3);
                R_z = T_z(1:3,1:3);
                
                sin_alpha = sin(v_i(4));
                cos_alpha = cos(v_i(4));
                sin_beta = sin(v_i(5));
                cos_beta = cos(v_i(5));
                sin_gamma = sin(v_i(6));
                cos_gamma = cos(v_i(6));
               
                R11 = R_i(1, 1);
                R12 = R_i(1, 2);
                R13 = R_i(1, 3);
                R21 = R_i(2, 1);
                R22 = R_i(2, 2);
                R23 = R_i(2, 3);
                R31 = R_i(3, 1);
                R32 = R_i(3, 2);
                R33 = R_i(3, 3);
                
                dR_alpha = [0 R13 -R12
                            0 R23 -R22
                            0 R33 -R32];
                dR_beta = [-sin_beta * cos_gamma   -sin_alpha * cos_beta * cos_gamma   -cos_alpha * cos_beta * cos_gamma
                           -sin_beta * sin_gamma   -sin_alpha * cos_beta * sin_gamma   -cos_alpha * cos_beta * sin_gamma
                               cos_beta                  -sin_alpha * sin_beta             -cos_alpha * sin_beta        ];
                dR_gamma = [-R21 -R22 -R23
                             R11  R12  R13
                              0    0    0 ];
                
                dt_ij = T_j(1:3, 4) - T_i(1:3, 4);
                
                R_z_R_i_ = [R_z'*dR_alpha'*dt_ij  R_z'*dR_beta'*dt_ij  R_z'*dR_gamma'*dt_ij];
                % Caluclate jacobians
                A = [-R_z'*R_i' R_z_R_i_; 
                      zeros(3)   -eye(3)];
                B = [R_z'*R_i' zeros(3); 
                     zeros(3)   eye(3) ];
                
                Jacobian(j_ind, i_ind) = A;
                Jacobian(j_ind, j_ind) = B;
                
                % Calculate error vector
                e = m2v((T_z) \ (T_i \ T_j)); 
                cost = cost + e'* omega * e;
                f_(j_ind) = e;
                % Formulate blocks
                H_ii =  A' * omega * A;
                H_ij =  A' * omega * B;
                H_jj =  B' * omega * B;
                b_i  = -A' * omega * e;
                b_j  = -B' * omega * e;
                
                % Update H and b matrix
                obj.H(i_ind,i_ind) = obj.H(i_ind,i_ind) + H_ii;
                obj.H(i_ind,j_ind) = obj.H(i_ind,j_ind) + H_ij;
                obj.H(j_ind,i_ind) = obj.H(j_ind,i_ind) + H_ij';
                obj.H(j_ind,j_ind) = obj.H(j_ind,j_ind) + H_jj;
                obj.b(i_ind) = obj.b(i_ind) + b_i;
                obj.b(j_ind) = obj.b(j_ind) + b_j;
            end
            obj.J = Jacobian;
      %      fprintf('det(J) = %d\n', det(Jacobian));
            obj.E = cost;
            obj.f = f_;
            fprintf('The loop detection cost = %d\n', cost);

        end

        function solve(obj)
            % Solves the linear system and update all pose node
         %   fprintf('Pose: %d, Edge: %d\n', obj.n_node, obj.n_edge);
            % The system (H b) is obtained only from relative constraints.
            % H is not full rank.
            % We solve this by anchoring the position of the 1st vertex
            % This can be expressed by adding teh equation
            % dx(1:3,1) = 0
            % which is equivalent to the following
            obj.H(1:6,1:6) = obj.H(1:6,1:6) + eye(6);
            det(obj.H);
            obj.H = obj.H + (obj.lambda) * eye(size(obj.H));
       %     H_sparse = sparse(obj.H);
            
            dx = pinv(obj.H) * obj.b;
  %          dx = inv(H_sparse + (obj.lambda) * eye(size(H_sparse))) * obj.b;
   %         dx = (H_sparse + (obj.lambda) * eye(size(H_sparse))) \ obj.b;
   %         plot(dx , '*-');
            obj.q = obj.J * dx;
            dpose = reshape(dx, 6, obj.n_node);
            
            % Update node with solution
            for i_node = 1:obj.n_node
                obj.node(i_node).pose = obj.node(i_node).pose ...
                    + dpose(:,i_node);
            end
        end
        
        function n_node = get.n_node(obj)
            n_node = numel(obj.node);
        end
        
        function n_edge = get.n_edge(obj)
            n_edge = numel(obj.edge);
        end
        
        function pose = get.pose(obj)
            pose = [obj.node.pose];
        end
        
    end  % methods public
    
end  % classdef


function ind = id2ind(id)
%ID2IND converts id to indices in H and b
ind = (6*(id-1)+1):(6*id);
end

 