% Computer 3D point cloud by triangulation

function [Xw, Yw, Zw] = Compute3D(phase, Mc, Mp, col, row, flag)

if (1 == flag) % h direction for pattern    
    Xw = zeros(row, col);
    Yw = zeros(row, col);
    Zw = zeros(row, col);
    f = (2 * pi ) / row;
    for i = 1 : 1 : row        
        for j = 1 : 1 : col  
            if(phase(i, j) == 0)
                continue;
            end
            xc = j;
            yc = i;
            yp = phase(yc, xc) / f;
            
            if( 0 == yp )
                continue;
            end

            C(1, 1) = Mc(1, 1) - Mc(3, 1) * xc;
            C(1, 2) = Mc(1, 2) - Mc(3, 2) * xc;
            C(1, 3) = Mc(1, 3) - Mc(3, 3) * xc;
            C(2, 1) = Mc(2, 1) - Mc(3, 1) * yc;
            C(2, 2) = Mc(2, 2) - Mc(3, 2) * yc;
            C(2, 3) = Mc(2, 3) - Mc(3, 3) * yc;
            C(3, 1) = Mp(2, 1) - Mp(3, 1) * yp;
            C(3, 2) = Mp(2, 2) - Mp(3, 2) * yp;
            C(3, 3) = Mp(2, 3) - Mp(3, 3) * yp;
            
            D(1) = Mc(3, 4) * xc - Mc(1, 4);
            D(2) = Mc(3, 4) * yc - Mc(2, 4);
            D(3) = Mp(3, 4) * yp - Mp(2, 4);
            
            world = Func_Inverse3(C) * D';

            Xw(yc, xc) = world(1);
            Yw(yc, xc) = world(2);
            Zw(yc, xc) = world(3);
        end   
    end
end

if (3 == flag) % h direction for pattern    
    Xw = 0;
    Yw = 0;
    Zw = 0;
    f = (2 * pi ) / row; 
    k = 1;
    for i = 1 : 3 : row  
        
        for j = 1 : 3 : col 
            if(phase(i, j) == 0)
                continue;
            end
            xc = j;
            yc = i;
            yp = phase(yc, xc) / f;
            
            if( 0 == yp )
                continue;
            end

            C(1, 1) = Mc(1, 1) - Mc(3, 1) * xc;
            C(1, 2) = Mc(1, 2) - Mc(3, 2) * xc;
            C(1, 3) = Mc(1, 3) - Mc(3, 3) * xc;
            C(2, 1) = Mc(2, 1) - Mc(3, 1) * yc;
            C(2, 2) = Mc(2, 2) - Mc(3, 2) * yc;
            C(2, 3) = Mc(2, 3) - Mc(3, 3) * yc;
            C(3, 1) = Mp(2, 1) - Mp(3, 1) * yp;
            C(3, 2) = Mp(2, 2) - Mp(3, 2) * yp;
            C(3, 3) = Mp(2, 3) - Mp(3, 3) * yp;
            
            D(1) = Mc(3, 4) * xc - Mc(1, 4);
            D(2) = Mc(3, 4) * yc - Mc(2, 4);
            D(3) = Mp(3, 4) * yp - Mp(2, 4);
            
            world = Func_Inverse3(C) * D';   
            Xw(k) = world(1);
            Yw(k) = world(2);
            Zw(k) = world(3);
            k = k + 1;
        end   
    end
end
if (0 == flag)  % v deriction
    Xw = zeros(row, col);
    Yw = zeros(row, col);
    Zw = zeros(row, col);
    f = (2 * pi ) / col;   
    for i = 1 : 1 :col  
        for j = 1 : 1 : row  
            xc = i;
            yc = j;
            if(phase(j, i) == 0)
                continue;
            end
            phi = phase(yc, xc); 
            xp = (2 * pi - phi) / f;
            
            if( 0 == xp )
                continue;
            end

            C(1, 1) = Mc(1, 1) - Mc(3, 1) * xc;
            C(1, 2) = Mc(1, 2) - Mc(3, 2) * xc;
            C(1, 3) = Mc(1, 3) - Mc(3, 3) * xc;
            C(2, 1) = Mc(2, 1) - Mc(3, 1) * yc;
            C(2, 2) = Mc(2, 2) - Mc(3, 2) * yc;
            C(2, 3) = Mc(2, 3) - Mc(3, 3) * yc;
            C(3, 1) = Mp(1, 1) - Mp(3, 1) * xp;
            C(3, 2) = Mp(1, 2) - Mp(3, 2) * xp;
            C(3, 3) = Mp(1, 3) - Mp(3, 3) * xp;
            
            D(1) = Mc(3, 4) * xc - Mc(1, 4);
            D(2) = Mc(3, 4) * yc - Mc(2, 4);
            D(3) = Mp(3, 4) * xp - Mp(1, 4);
            
            world = Func_Inverse3(C) * D';

            Xw(yc, xc) = world(1);
            Yw(yc, xc) = world(2);
            Zw(yc, xc) = world(3);
        end
    end
end

if (2 == flag)  % v deriction pattern
    f = (2 * pi ) / col;
    k = 1;
    for i = 1 : 6 :col         
        for j = 1 : 6 : row  
            xc = i;
            yc = j;
            if(phase(j, i) == 0)
                continue;
            end
            phi = phase(yc, xc); 
            xp = (2 * pi - phi) / f;
            
            if( 0 == xp )
                continue;
            end

            C(1, 1) = Mc(1, 1) - Mc(3, 1) * xc;
            C(1, 2) = Mc(1, 2) - Mc(3, 2) * xc;
            C(1, 3) = Mc(1, 3) - Mc(3, 3) * xc;
            C(2, 1) = Mc(2, 1) - Mc(3, 1) * yc;
            C(2, 2) = Mc(2, 2) - Mc(3, 2) * yc;
            C(2, 3) = Mc(2, 3) - Mc(3, 3) * yc;
            C(3, 1) = Mp(1, 1) - Mp(3, 1) * xp;
            C(3, 2) = Mp(1, 2) - Mp(3, 2) * xp;
            C(3, 3) = Mp(1, 3) - Mp(3, 3) * xp;
            
            D(1) = Mc(3, 4) * xc - Mc(1, 4);
            D(2) = Mc(3, 4) * yc - Mc(2, 4);
            D(3) = Mp(3, 4) * xp - Mp(1, 4);
            
            world = Func_Inverse3(C) * D';

            Xw(k) = world(1);
            Yw(k) = world(2);
            Zw(k) = world(3);
            k = k + 1;
        end
    end
end

end
