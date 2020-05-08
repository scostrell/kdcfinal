function N = dV_barrett(th1, th2, th3, th4, th5, th6, th7, m1, m2, m3, m4, m5, m6, m7, gth_sl1, gth_sl2, gth_sl3,gth_sl4,gth_sl5,gth_sl6,gth_sl7)
% A function to compute the derivative of potential energy due to gravity
% See pp. 169, 174 of MLS
    g = 9.81;
    n = 7;
    
    thetas = [th1, th2, th3, th4, th5, th6, th7];
    masses = [m1, m2, m3, m4, m5, m6, m7];
    h = [gth_sl1(3,4), gth_sl2(3,4), gth_sl3(3,4),gth_sl4(3,4),...
        gth_sl5(3,4),gth_sl6(3,4),gth_sl7(3,4)];
    V = 0;
    for i = 1:n
        V = V + g*masses(i)*h(i);
    end
    
    N = [];
    for i = 1:n
        N = [N;diff(V,thetas(i))];
    end

end