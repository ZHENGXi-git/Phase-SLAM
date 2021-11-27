function R_t = v2m_v(v)

c = cos(v(4));
s = sin(v(4));
% 注意 R_T的格式！

R_t = [1   0   0   0
       0   c   -s  v(2)
       0   s   c  v(3)
       0   0   0   1];
end