function iA = inv_rank1_add( iA, x )


y = iA*x;

%iA = iA - (1/(1+x'*y))* (y* x'*iA); 
%iA = iA - (1/(1+x'*y))* (y* y'); 

iA = iA -  ( (1/(1+x'*y))* y) * y';