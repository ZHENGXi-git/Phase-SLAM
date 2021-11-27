%--------------point cloud fusion--------------%
function[X_, Y_, Z_] = PCFusion(X, Y, Z, R_t, times, flag)

[row, col] = size(X);
X_ = zeros(row, col);
Y_ = zeros(row, col);
Z_ = zeros(row, col);
width = col;
height = row;

if (1 == flag)
    R = R_t(1:3, 1:3);
    t = R_t(1:3, 4);
    for i = 1: 1 : height
        for j = 1: 1: width
            xc = j;
            yc = i;
            x = X(yc, xc);  
            y = Y(yc, xc);
            z = Z(yc, xc);
            Pos = R * [x y z]' + t;
            X(yc, xc) = Pos(1);
            Y(yc, xc) = Pos(2);
            Z(yc, xc) = Pos(3);
        end
    end
    X_ = X;
    Y_ = Y;
    Z_ = Z;
end

if (2 == flag)
    for k = times: -1: 1
        R = R_t(k).R;
        T = R_t(k).t;
        for i = 1: 1 : height
            for j = 1: 1: width
                xc = j;
                yc = i;
                x = X(yc, xc); 
                y = Y(yc, xc);
                z = Z(yc, xc);                
                Pos = R' * ([x y z]' - T);               
                X(yc, xc) = Pos(1);
                Y(yc, xc) = Pos(2);
                Z(yc, xc) = Pos(3);
            end
        end
    end
    X_ = X;
    Y_ = Y;
    Z_ = Z;
end

if (3 == flag)
    R = R_t(1:3, 1:3);
    t = R_t(1:3, 4);
    for i = 1: 1 : height
        for j = 1: 1: width
            xc = j;
            yc = i;
            x = X(yc, xc); 
            y = Y(yc, xc);
            z = Z(yc, xc);
            Pos = R' * ([x y z]' - t);
            X(yc, xc) = Pos(1);
            Y(yc, xc) = Pos(2);
            Z(yc, xc) = Pos(3);
        end
    end
    X_ = X;
    Y_ = Y;
    Z_ = Z;
end

if (4 == flag)
    for k = 1: 1: times
        R = R_t(k).R;
        t = R_t(k).t;
        for i = 1: 1 : height
            for j = 1: 1: width
                xc = j;
                yc = i;
                x = X(yc, xc); 
                y = Y(yc, xc);
                z = Z(yc, xc);
                Pos = R * [x y z]' + t;
                X(yc, xc) = Pos(1);
                Y(yc, xc) = Pos(2);
                Z(yc, xc) = Pos(3);
            end
        end
    end
    X_ = X;
    Y_ = Y;
    Z_ = Z;
end

end