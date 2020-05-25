function [cor overlap ind  SSRs] = c_corr( C1, C2, Thresh,   G, D1, D2) 

% computes the maximum correlations between the columns of C1 and C2


if nargin < 3 || isempty(Thresh),  Thresh = 1e-3; end


[m n] = size(C1); if m<n, C1 = C1'; end
 
[m n] = size(C2); if m<n, C2 = C2'; end

if ~isequal(size(C1),size(C1)), 
    fprintf( 'Different Sizes\n' );
    return
end 

ind = zeros( size(C1,2), 1 ); 
cor = zeros( size(ind) );
overlap = zeros( size(ind,1), 3 );

r = corr( C1, C2 );
for i=1:size(r,2)
    [cor(i) ind(i)] = max( abs(r(:,i)) );
    
    if Thresh < 1
        c1 = double( abs( C1(:,ind(i)) ) > Thresh ); 
        c2 = double( abs( C2(:,i) ) > Thresh ); 
    else 
        vals = sort( abs(C2(:,i)) ); 
        thresh = vals(end-Thresh);
        c1 = double( abs( C1(:,ind(i)) ) > 0 );
        c2 = double( abs( C2(:,i) ) > thresh );
    end 
    
    overlap(i,1) = sum( c1 );
    overlap(i,2) = sum( c2 );
    overlap(i,3) = c1' * c2;
end 


 
% Computes the sum of squared residuals

if nargin >=6 && nargout >=4 && ~isempty(G)
 SSRs = zeros( 4, 1 );
 SSRs(1) = norm( G, 'fro' );
 SSRs(2) = norm( G - D1*C1',    'fro' ); % Data with (G) and without noise (D*C)  
 SSRs(3) = norm( G - D2*C2',    'fro' ); 
 SSRs(4) = norm( D1*C1'-D2*C2', 'fro' ); % Data without noise and model prediction  
end  











