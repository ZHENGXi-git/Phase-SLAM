function R_t = v2m(v)

R_t = zeros(4);
cos_alpha = cos(v(4));
sin_alpha = sin(v(4));
cos_beta = cos(v(5));
sin_beta = sin(v(5));
cos_gamma = cos(v(6));
sin_gamma = sin(v(6));

R11 = cos_beta * cos_gamma;
R12 = -sin_alpha * sin_beta * cos_gamma - cos_alpha * sin_gamma;
R13 = -cos_alpha * sin_beta * cos_gamma + sin_alpha * sin_gamma;
R21 = cos_beta * sin_gamma;
R22 = -sin_alpha * sin_beta * sin_gamma + cos_alpha * cos_gamma;
R23 = -cos_alpha * sin_beta * sin_gamma - sin_alpha * cos_gamma;
R31 = sin_beta;
R32 = sin_alpha * cos_beta;
R33 = cos_alpha * cos_beta;

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
t = [v(1); v(2); v(3)];
h = [0 0 0 1];
R_t = [R, t];
R_t = [R_t; h];

% 注意 R_T的格式！
R_t = [ R11  R12  R13  v(1)  
        R21  R22  R23  v(2)  
        R31  R32  R33  v(3)
         0    0    0    1  ];
     
end