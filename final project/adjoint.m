function ad = adjoint(g)
    R = g(1:3,1:3);
    p = g(1:3,4);
    phat = matcross(p);
    ad = [R phat*R; zeros(3) R];
end