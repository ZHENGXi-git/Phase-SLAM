function[pose, Rn, tn, R_t] = PoseCompution(motion, Rn, tn, R_t, count, flag)

cos_alpha = cos(motion(4)); 
sin_alpha = sin(motion(4));
cos_beta = cos(motion(5));
sin_beta = sin(motion(5));
cos_gamma = cos(motion(6));
sin_gamma = sin(motion(6));

R1 = [1        0      0
    0 cos_alpha -sin_alpha
    0 sin_alpha cos_alpha];
R2 = [cos_beta 0 -sin_beta;
    0    1    0;
    sin_beta 0 cos_beta];
R3 = [cos_gamma -sin_gamma 0;
    sin_gamma cos_gamma 0;
    0        0      1];
R = R3 * R2 * R1;
t = [motion(1); motion(2); motion(3)];
R_t(count).R = R;
R_t(count).t = t;
Rn = Rn * R';
tn = -Rn * t + tn ;
[alpha, beta, gamma] = dcm2angle(Rn, 'XYZ');
if (0 == flag)
    % position
    [px, py, pz] = PCFusion(0, 0, 0, R_t, count, 2);
    % rotation
    pose = [px, py, pz, -alpha, beta, -gamma];
end
if (1 == flag)
    [px, py, pz] = PCFusion(0, 0, 0, R_t, count, 4);
    pose = [py, pz, alpha];
end

end