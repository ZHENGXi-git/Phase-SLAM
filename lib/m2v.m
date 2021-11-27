function vec = m2v(R_t)
vec(1, 1) = R_t(1, 4);
vec(2, 1) = R_t(2, 4);
vec(3, 1) = R_t(3, 4);
R = R_t(1: 3, 1: 3);
[alpha, beta, gamma] = dcm2angle(R, 'XYZ'); 
vec(4, 1) = -alpha;
vec(5, 1) = beta;
vec(6, 1) = -gamma;
end