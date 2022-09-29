function print_array(M)
fprintf('double M_expert[N_expert][4][4] = {{{');
for i=1:4
    for j=1:4
        fprintf('%.16f', M(i,j));
        if j<4, fprintf(', '); end
    end
    if i<4, fprintf('},\n                                    {');
    else, fprintf('}}};\n'); end
end
end