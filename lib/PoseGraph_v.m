classdef PoseGraph_v < handle
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
        
        function obj = PoseGraph_v()
            % Constructor of PoseGraph
            obj.node = PoseNode_v.empty;
            obj.edge = PoseEdge_v.empty;
        end
        
        function readGraph_v(obj, vertices, edges)
                  
          %  vertices = vertice.nodes;
            for i_node = 1:size(vertices,2)
                vi = vertices(:,i_node);
                id = vi(1) + 1;
                pose_(1) = 0;      %x
                pose_(2) = vi(2);  %y
                pose_(3) = vi(3);  %z
                pose_(4) = vi(4);  %alpha
                
                obj.node(i_node) = PoseNode_v(id, pose_);  % y z alpha
            end

       %     edges = edge.edge;
            for i_edge = 1:size(edges,2)
                ei = edges(:,i_edge);
                id_from   = ei(1) + 1;
                id_to     = ei(2) + 1;
                mean(1)   = 0;
                mean(2)   = ei(3);
                mean(3)   = ei(4);
                mean(4)   = ei(5);

                infm = eye(4,4);
            %    infm
                obj.edge(i_edge) = PoseEdge_v(id_from, id_to, mean, infm);
            end

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
            obj.H = zeros(obj.n_node*4);   % 3n x 3n square matrix
            obj.b = zeros(obj.n_node*4,1); % 3n x 1  column vector
            
      %      fprintf('Linearizing.\n');
            obj.linearize_();
            
      %      fprintf('Solving.\n');
            obj.solve();
        end

        function linearize_(obj)

            f_ = zeros(obj.n_node*4, 1); 
            cost = 0;
            obj.E = 0;
            Jacobian = zeros(obj.n_node*4, obj.n_node*4);
            for i_edge = 1:obj.n_edge
            
                ei = obj.edge(i_edge);   %edge
                % Get edge information
                j_node = ei.id_from;
                i_node = ei.id_to;
                T_z = v2m_v(ei.mean);
                omega = ei.infm; 
                
                % Get node information
                v_i = obj.node(i_node).pose;
                v_j = obj.node(j_node).pose;
                
                i_ind = id2ind(i_node);    % H_spares
                j_ind = id2ind(j_node);
                
                T_i = v2m_v(v_i);
                T_j = v2m_v(v_j);
                R_i = T_i(1:3,1:3);
                R_z = T_z(1:3,1:3);
                
                s = sin(v_i(4));
                c = cos(v_i(4));
              
                dR_i = [1  0  0
                        0 -s -c
                        0  c -s];
                dt_ij = T_j(1:3, 4) - T_i(1:3, 4);
                
                % Caluclate jacobians
                A = [-R_z'*R_i' R_z'*dR_i'*dt_ij; 0 0 0 -1];
                B = [R_z'*R_i' [0 0 0]'; 0 0 0 1];
                
                Jacobian(j_ind, i_ind) = A;
                Jacobian(j_ind, j_ind) = B;
                
                % Calculate error vector
                e = m2v_v(inv(T_z) * inv(T_i) * T_j);  
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
            obj.H(1:4,1:4) = obj.H(1:4,1:4) + eye(4);
            det(obj.H);
            obj.H = obj.H + (obj.lambda) * eye(size(obj.H));
       %     H_sparse = sparse(obj.H);
            
            dx = pinv(obj.H) * obj.b;
  %          dx = inv(H_sparse + (obj.lambda) * eye(size(H_sparse))) * obj.b;
   %         dx = (H_sparse + (obj.lambda) * eye(size(H_sparse))) \ obj.b;
   %         plot(dx , '*-');
            obj.q = obj.J * dx;
            dpose = reshape(dx, 4, obj.n_node);
            
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
ind = (4*(id-1)+1):(4*id);
end

 