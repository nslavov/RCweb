function iA = inv_rank1( iA, x, TYPE )



switch TYPE
    case { 1, 'red' }

    y = iA*x;

    iA = iA + (1/(1-x'*y))* (y* x'*iA);
    
    case { 2, 'add' }
        
    y = iA*x;

    iA = iA - (1/(1+x'*y))* (y* x'*iA);
end 