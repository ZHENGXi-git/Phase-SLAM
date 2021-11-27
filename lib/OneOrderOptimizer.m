%-------------------------------------%
%Pose Estimation Optimizer
%-------------------------------------%

function[motion, diff] = OneOrderOptimizer(phase1, phase2, X, Y, Z, ROI, Kp, Mc, initialization, iteration, flag)

initialization(4: 6) = initialization(4: 6) / 180 * pi;

motion = initialization(1: 6);
state = initialization(1: 6);  
diff = 10000000;
[row, col] = size(phase1);

grad_v = zeros(row, col);
grad_u = zeros(row, col);
cost = zeros(row, col);

point_col = ROI(1); %x
point_row = ROI(2); %y
height = ROI(3);    
width = ROI(4);   

eta = initialization(7: 12);

for i = 2 : 1 : row - 1
    dy_1 = phase2(i + 1, :) - phase2(i, :);  
    dy_2 = phase2(i, :) - phase2(i - 1, :);
    dy_3 = phase2(i + 1, :) - phase2(i - 1, :);
    grad_v(i,:) = dy_1 / 4 + dy_2 / 4 + dy_3 / 2;
end
for j = 2 : 1 : col - 1
    dx_1 = phase2(:, j + 1) - phase2(:, j);  
    dx_2 = phase2(:, j) - phase2(:, j - 1);
    dx_3 = phase2(:, j + 1) - phase2(:, j - 1);
    grad_u(:, j) = dx_1 / 4 + dx_2 / 4 + dx_3 / 2;
end


if (0 == flag)  % horizontal

    K = row / (2 * pi);
    fyp = Kp(2, 2);
    Cyp = Kp(2, 3);
    Phase = zeros(row, col);
    
    for k = 1: 1: iteration

        J = zeros(6, 1);
        g = zeros(6, 1);
        delta = zeros(6, 1);
        
        COST = 0;
        cos_alpha = cos(state(4));
        sin_alpha = sin(state(4));
        cos_beta = cos(state(5));
        sin_beta = sin(state(5));
        cos_gamma = cos(state(6));
        sin_gamma = sin(state(6));
        R11 = cos_beta * cos_gamma;
        R12 = -sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma;
        R13 = -cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma;
        R21 = cos_beta * sin_gamma;
        R22 = -sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma;
        R23 = -cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma;
        R31 = sin_beta;
        R32 = sin_alpha * cos_beta;
        R33 = cos_alpha * cos_beta;
        
        for i = 1: 1: height
            v = point_row + i;
            for j = 1: 1: width
                u = point_col + j;
                z = Z(v, u);
                y = Y(v, u);
                x = X(v, u);
                if (abs(z) < 0.001)
                    continue;
                end

                x_ = R11 * x + R12 * y + R13 * z + state(1);
                y_ = R21 * x + R22 * y + R23 * z + state(2);
                z_ = R31 * x + R32 * y + R33 * z + state(3);
                s_ =  Mc(3, 1) * x_ + Mc(3, 2) * y_ + Mc(3, 3) * z_ + Mc(3, 4);
                u_ = (Mc(1, 1) * x_ + Mc(1, 2) * y_ + Mc(1, 3) * z_ + Mc(1, 4)) / s_;  % xc_
                v_ = (Mc(2, 1) * x_ + Mc(2, 2) * y_ + Mc(2, 3) * z_ + Mc(2, 4)) / s_;  % yc_
                
                if(u_ < 1 || v_ < 1 || u_ > (col-1) || v_ > (row-1))
                    continue;
                end
                x0 = floor(u_);
                x1 = x0 + 1;
                y0 = floor(v_);
                y1 = y0 + 1;
                weight_x = u_ - x0;
                weight_y = v_ - y0;

                phi_2 = phase2(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + phase2(y1, x0) * weight_y * (1 - weight_x)...
                    + phase2(y0, x1) * (1 - weight_y) * weight_x...
                    + phase2(y1, x1) * weight_x * weight_y;
                gradu = grad_u(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + grad_u(y1, x0) * weight_y * (1 - weight_x)...
                    + grad_u(y0, x1) * (1 - weight_y) * weight_x...
                    + grad_u(y1, x1) * weight_x * weight_y;
                gradv = grad_v(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + grad_v(y1, x0) * weight_y * (1 - weight_x)...
                    + grad_v(y0, x1) * (1 - weight_y) * weight_x...
                    + grad_v(y1, x1) * weight_x * weight_y;
                
                phi_ = fyp * y_ / (z_ * K) + Cyp / K;
                Phase(y0, x0) = phi_;
                cost(v, u) = phi_ - phi_2;
                COST = COST + 0.5 * cost(v, u) * cost(v, u);
                
                J_x_x = 1;
                J_x_y = 0;
                J_x_z = 0;
                J_y_x = 0;
                J_y_y = 1;
                J_y_z = 0;
                J_z_x = 0;
                J_z_y = 0;
                J_z_z = 1;
                J_x_alpha = R13 * y - R12 * z;
                J_x_beta  = -sin_beta * cos_gamma * x + (-sin_alpha * cos_beta * cos_gamma) * y + (-cos_alpha * cos_beta * cos_gamma) * z;
                J_x_gamma = state(2) - y_;
                J_y_alpha = R23 * y - R22 * z;
                J_y_beta  = -sin_beta * sin_gamma * x + (-sin_alpha * cos_beta * sin_gamma) * y + (-cos_alpha * cos_beta * sin_gamma) * z;
                J_y_gamma = x_ - state(1);
                J_z_alpha = R33 * y - R32 * z;
                J_z_beta  = cos_beta * x + (-sin_alpha * sin_beta) * y + (-cos_alpha * sin_beta) * z;
                J_z_gamma = 0;
                a_u = Mc(1, 1) - Mc(3, 1) * u_;
                b_u = Mc(1, 2) - Mc(3, 2) * u_;
                c_u = Mc(1, 3) - Mc(3, 3) * u_;
                J_u_x = (a_u * J_x_x + b_u * J_y_x + c_u * J_z_x) / s_; % = a_u / s_
                J_u_y = (a_u * J_x_y + b_u * J_y_y + c_u * J_z_y) / s_; % = b_u / s_
                J_u_z = (a_u * J_x_z + b_u * J_y_z + c_u * J_z_z) / s_;
                J_u_alpha = (a_u * J_x_alpha + b_u * J_y_alpha + c_u * J_z_alpha) / s_;
                J_u_beta  = (a_u * J_x_beta  + b_u * J_y_beta  + c_u * J_z_beta) / s_;
                J_u_gamma = (a_u * J_x_gamma + b_u * J_y_gamma + c_u * J_z_gamma) / s_;
                a_v = Mc(2, 1) - Mc(3, 1) * v_;
                b_v = Mc(2, 2) - Mc(3, 2) * v_;
                c_v = Mc(2, 3) - Mc(3, 3) * v_;
                J_v_x = (a_v * J_x_x + b_v * J_y_x + c_v * J_z_x) / s_; % = a_v / s_
                J_v_y = (a_v * J_x_y + b_v * J_y_y + c_v * J_z_y) / s_; % b_v / s_
                J_v_z = (a_v * J_x_z + b_v * J_y_z + c_v * J_z_z) / s_; % c_v / s_
                J_v_alpha = (a_v * J_x_alpha + b_v * J_y_alpha + c_v * J_z_alpha) / s_;
                J_v_beta  = (a_v * J_x_beta  + b_v * J_y_beta  + c_v * J_z_beta) / s_;
                J_v_gamma = (a_v * J_x_gamma + b_v * J_y_gamma + c_v * J_z_gamma) / s_;
                
                J_phi_x = fyp / K * (J_y_x * z_ - J_z_x * y_) / z_^2;   % = 0
                J_phi_y = fyp / K * (J_y_y * z_ - J_z_y * y_) / z_^2;   % J_z_y = 0
                J_phi_z = fyp / K * (J_y_z * z_ - J_z_z * y_) / z_^2;   % J_y_z = 0
                J_phi_alpha = fyp / K * (J_y_alpha * z_ - J_z_alpha * y_) / z_^2;
                J_phi_beta  = fyp / K * (J_y_beta * z_  - J_z_beta * y_) / z_^2;
                J_phi_gamma = fyp / K * (J_y_gamma * z_ - J_z_gamma * y_) / z_^2;

                J(1) = J_phi_x - (gradu * J_u_x + gradv * J_v_x);
                J(2) = J_phi_y - (gradu * J_u_y + gradv * J_v_y);
                J(3) = J_phi_z - (gradu * J_u_z + gradv * J_v_z);
                J(4) = J_phi_alpha - (gradu * J_u_alpha + gradv * J_v_alpha);
                J(5) = J_phi_beta  - (gradu * J_u_beta  + gradv * J_v_beta);
                J(6) = J_phi_gamma - (gradu * J_u_gamma + gradv * J_v_gamma);
                         
                g = g + cost(v, u) * J;
                
            end
        end        
        diff = COST;         
        delta = -eta .* g';

        if(abs(delta(1)) < 0.0005 && abs(delta(2)) < 0.0005 && abs(delta(3)) < 0.0005 && abs(delta(4) * 180 / pi) < 0.0001 && abs(delta(5) * 180 / pi) < 0.0001 && abs(delta(6) * 180 / pi) < 0.0001)  %% 多变量模式 是否要考虑当其中一个变量步长很小时，是否需要锁定？
            motion = state; 
            break;
        end
        state = state + delta;
        motion = state;       
    end
    
end

if (1 == flag) % vertical
    K = col / (2 * pi);
    fxp = Kp(1, 1);
    alpha_c = Kp(1, 2);
    Cxp = Kp(1, 3);

    for k = 1: 1: iteration
        
        J = zeros(6, 1);
        g = zeros(6, 1);
        delta = zeros(6, 1);
        
        COST = 0;
        cos_alpha = cos(state(4));
        sin_alpha = sin(state(4));
        cos_beta = cos(state(5));
        sin_beta = sin(state(5));
        cos_gamma = cos(state(6));
        sin_gamma = sin(state(6));
        R11 = cos_beta * cos_gamma;
        R12 = -sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma;
        R13 = -cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma;
        R21 = cos_beta * sin_gamma;
        R22 = -sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma;
        R23 = -cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma;
        R31 = sin_beta;
        R32 = sin_alpha * cos_beta;
        R33 = cos_alpha * cos_beta;
        
        for i = 1: 1: height
            v = point_row + i;
            for j = 1: 1: width
                u = point_col + j;
                z = Z(v, u);
                y = Y(v, u);
                x = X(v, u);
                if (abs(z) < 0.001)
                    continue;
                end
                x_ = R11 * x + R12 * y + R13 * z + state(1);
                y_ = R21 * x + R22 * y + R23 * z + state(2);
                z_ = R31 * x + R32 * y + R33 * z + state(3);
                s_ =  Mc(3, 1) * x_ + Mc(3, 2) * y_ + Mc(3, 3) * z_ + Mc(3, 4);
                u_ = (Mc(1, 1) * x_ + Mc(1, 2) * y_ + Mc(1, 3) * z_ + Mc(1, 4)) / s_;  % xc_
                v_ = (Mc(2, 1) * x_ + Mc(2, 2) * y_ + Mc(2, 3) * z_ + Mc(2, 4)) / s_;  % yc_
                
                if(u_ < 1 || v_ < 1 || u_ > (col-1) || v_ > (row-1))
                    continue;
                end
                x0 = floor(u_);
                x1 = x0 + 1;
                y0 = floor(v_);
                y1 = y0 + 1;
                weight_x = u_ - x0;
                weight_y = v_ - y0;
                phi_2 = phase2(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + phase2(y1, x0) * weight_y * (1 - weight_x)...
                    + phase2(y0, x1) * (1 - weight_y) * weight_x...
                    + phase2(y1, x1) * weight_x * weight_y;
                gradu = grad_u(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + grad_u(y1, x0) * weight_y * (1 - weight_x)...
                    + grad_u(y0, x1) * (1 - weight_y) * weight_x...
                    + grad_u(y1, x1) * weight_x * weight_y;
                gradv = grad_v(y0, x0) * (1 - weight_x) * (1 - weight_y)...
                    + grad_v(y1, x0) * weight_y * (1 - weight_x)...
                    + grad_v(y0, x1) * (1 - weight_y) * weight_x...
                    + grad_v(y1, x1) * weight_x * weight_y;

                phi_ = 2 * pi - ((fxp * x_ + alpha_c * y_) / (K * z_) + Cxp / K);
                cost(v, u) = phi_ - phi_2;
                COST = COST + 0.5 * cost(v, u) * cost(v, u);
                
                J_x_x = 1;
                J_x_y = 0;
                J_x_z = 0;
                J_y_x = 0;
                J_y_y = 1;
                J_y_z = 0;
                J_z_x = 0;
                J_z_y = 0;
                J_z_z = 1;
                J_x_alpha = R13 * y - R12 * z;
                J_x_beta  = -sin_beta * cos_gamma * x + (-sin_alpha * cos_beta * cos_gamma) * y + (-cos_alpha * cos_beta * cos_gamma) * z;
                J_x_gamma = -y_;
                J_y_alpha = R23 * y - R22 * z;
                J_y_beta  = -sin_beta * sin_gamma * x + (-sin_alpha * cos_beta * sin_gamma) * y + (-cos_alpha * cos_beta * sin_gamma) * z;
                J_y_gamma = x_;
                J_z_alpha = R33 * y - R32 * z;
                J_z_beta  = cos_beta * x + (-sin_alpha * sin_beta) * y + (-cos_alpha * sin_beta) * z;
                J_z_gamma = 0;
                a_u = Mc(1, 1) - Mc(3, 1) * u_;
                b_u = Mc(1, 2) - Mc(3, 2) * u_;
                c_u = Mc(1, 3) - Mc(3, 3) * u_;
                J_u_x = (a_u * J_x_x + b_u * J_y_x + c_u * J_z_x) / s_;
                J_u_y = (a_u * J_x_y + b_u * J_y_y + c_u * J_z_y) / s_;
                J_u_z = (a_u * J_x_z + b_u * J_y_z + c_u * J_z_z) / s_;
                J_u_alpha = (a_u * J_x_alpha + b_u * J_y_alpha + c_u * J_z_alpha) / s_;
                J_u_beta  = (a_u * J_x_beta  + b_u * J_y_beta  + c_u * J_z_beta) / s_;
                J_u_gamma = (a_u * J_x_gamma + b_u * J_y_gamma + c_u * J_z_gamma) / s_;
                a_v = Mc(2, 1) - Mc(3, 1) * v_;
                b_v = Mc(2, 2) - Mc(3, 2) * v_;
                c_v = Mc(2, 3) - Mc(3, 3) * v_;
                J_v_x = (a_v * J_x_x + b_v * J_y_x + c_v * J_z_x) / s_;
                J_v_y = (a_v * J_x_y + b_v * J_y_y + c_v * J_z_y) / s_;
                J_v_z = (a_v * J_x_z + b_v * J_y_z + c_v * J_z_z) / s_;
                J_v_alpha = (a_v * J_x_alpha + b_v * J_y_alpha + c_v * J_z_alpha) / s_;
                J_v_beta  = (a_v * J_x_beta  + b_v * J_y_beta  + c_v * J_z_beta) / s_;
                J_v_gamma = (a_v * J_x_gamma + b_v * J_y_gamma + c_v * J_z_gamma) / s_;
                
                J_phi_x = -((fxp * J_x_x + alpha_c * J_y_x) * z_ - J_z_x * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                J_phi_y = -((fxp * J_x_y + alpha_c * J_y_y) * z_ - J_z_y * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                J_phi_z = -((fxp * J_x_z + alpha_c * J_y_z) * z_ - J_z_z * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                J_phi_alpha = -((fxp * J_x_alpha + alpha_c * J_y_alpha) * z_ - J_z_alpha * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                J_phi_beta  = -((fxp * J_x_beta  + alpha_c * J_y_beta)  * z_ - J_z_beta  * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                J_phi_gamma = -((fxp * J_x_gamma + alpha_c * J_y_gamma) * z_ - J_z_gamma * (fxp * x_ + alpha_c * y_)) / (z_^2 * K);
                
                J(1) = J_phi_x - (gradu * J_u_x + gradv * J_v_x);
                J(2) = J_phi_y - (gradu * J_u_y + gradv * J_v_y);
                J(3) = J_phi_z - (gradu * J_u_z + gradv * J_v_z);
                J(4) = J_phi_alpha - (gradu * J_u_alpha + gradv * J_v_alpha);
                J(5)  = J_phi_beta  - (gradu * J_u_beta  + gradv * J_v_beta);
                J(6) = J_phi_gamma - (gradu * J_u_gamma + gradv * J_v_gamma);
                         
                g = g + cost(v, u) * J;
                
            end
        end
        diff = COST;
        delta = -eta .* g';
        if(abs(delta(1)) < 0.0001 && abs(delta(2)) < 0.0001 && abs(delta(3)) < 0.0001 && abs(delta(4) * 180 / pi) < 0.0001 && abs(delta(5) * 180 / pi) < 0.0001 && abs(delta(6) * 180 / pi) < 0.0001)  %% 多变量模式 是否要考虑当其中一个变量步长很小时，是否需要锁定？
            motion = state;  
            break;
        end
        state = state + delta;
        motion = state;
    end
    
end

end