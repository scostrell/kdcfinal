function N = dV(th1, th2, th3, m1, m2, m3, gth_sl1, gth_sl2, gth_sl3, g)
% A function to compute the derivative of potential energy due to gravity
% See pp. 169, 174 of MLS

    h1 = gth_sl1(3,4);
    h2 = gth_sl2(3,4);
    h3 = gth_sl3(3,4);
    
    V = g*(m1*h1 + m2*h2 + m3*h3);
    N = [diff(V,th1); diff(V,th2); diff(V,th3)];
end