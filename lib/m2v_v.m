function vec = m2v_v(R_t)
vec(1, 1) = 0;
vec(2, 1) = R_t(2, 4);
vec(3, 1) = R_t(3, 4);
vec(4, 1) = atan2(R_t(3, 2), R_t(2, 2));
end