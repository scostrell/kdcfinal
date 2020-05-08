function ad = adinv(g)
% p.56 MLS
    R = g(1:3,1:3);
    p = g(1:3,4);
    phat = matcross(p);
    ad = [R' -R' * phat; zeros(3) R'];
end