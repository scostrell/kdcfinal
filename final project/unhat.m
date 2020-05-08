function vec = unhat(mat)
    v1 = -mat(2,3);
    v2 = mat(1,3);
    v3 = -mat(1,2);
    
    vec = [v1 v2 v3];
%     vec_hat = [0 -v3 v2;
%                v3 0 -v1;
%                -v2 v1 0];
end