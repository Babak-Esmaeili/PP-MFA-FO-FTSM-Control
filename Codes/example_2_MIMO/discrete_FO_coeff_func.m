function out = discrete_FO_coeff_func(beta,L,h)

    out = zeros(1,L+1);
    out(1) = 1;
    
    vec = beta*ones(1,L) - (0:1:L-1);

    fact = 1;
    
    for j = 1:L
       
        fact = fact*j;
        out(j+1) = (h^(-beta))*(((-1)^j)*prod(vec(1:j)))/fact;
        
    end

end
