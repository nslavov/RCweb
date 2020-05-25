function [G V D C Non0 s time til_D B_real] = data_Gen( Param, ...
         TYPE, Num_Obs, N, P, Sparsness, Noise_in_C, Noise_in_G )
     
     
% Generates data from a sparse network that is used by SpAce to infer the
% net. data_Gen takes either a single structure of variables (Param) as an
% argument or multiple variables if Param is empty array.  
     


if ~isempty(Param)
    TYPE = Param.TYPE;
    Num_Obs = Param.Num_Obs;
    N = Param.N; 
    P = Param.P;
    Sparsness = Param.Sparsness;
    Noise_in_C = Param.Noise_in_C;
    Noise_in_G = Param.Noise_in_G;   
end 

D = rand( Num_Obs, P ); 

if  isempty(TYPE) || TYPE < 1, 
        C = rand( P, N );
        C = ( C > Sparsness ); 
        ind = find( sum(C)~=0 );
        for i=1:numel(ind)
           C( randi(P), ind(i) ) = 1; 
        end 
elseif TYPE > 1,
       [C sp deg] = gen_PowLaw_Net( N, P  );    
end 
C = (rand( size(C) )-0.5)  .* C;

Non0 = sum( abs(C)>0, 2 );

if  (~isempty(Param)|| nargin >= 5) && ~isempty(Noise_in_C) 
    C = C + Noise_in_C*( rand( size(C) ) - 0.5 ); 
end

G = D * C;

if  (~isempty(Param)|| nargin >= 6) && ~isempty(Noise_in_G) 
    Noise = (rand( size(G) )-0.5);
    G = G + ( Noise_in_G * var(G(:))/var(Noise(:)) ) * Noise;
end
 

tm = cputime;
[U S V] = svd( G, 'econ'  );  s = diag(S); 
time = cputime - tm; 

til_D = U * S;  % G * V(:,1:P);

if nargout >= 9
    a = inv( til_D' * til_D ) * (til_D' * D ); 

    if size(a,1)==size(a,2), 
        B_real = inv( a )';  
    end 
end 









function [C sp deg] = gen_PowLaw_Net( N, P )  %sparsity

% Sparsity must be less than 1  

C = rand( P, N );

%Thresh = linspace( 0, 1, P ); 
%Thresh = 2.^Thresh;

Thresh = logspace( -0.45, -0.04, P ); 

%Thresh = ( sparsity/0.43 ) * Thresh;

for i=1:P
    
    C(i,:) = C(i,:) > Thresh(i);  
    
end 
 
deg = sum( C, 2 ); 

sp = 1 - sum(deg)/(N*P);











function C = gen_Power_Law_Net( N, P, sparsity, Exp )


if nargin < 4 || isempty(Exp), Exp = 2;  end 

mn = min(P,N);
num = P*N; 

C = zeros( P, N ); 
seed = randperm( num ); 
C ( seed(1:N) )  = rand( N, 1 );  
Non_Zero  = round( (1-sparsity) * num ); 
k=N;
%%
Sum_2 = sum( C, 2 )+(N/Exp);
while k<Non_Zero

    p =  Sum_2  .* rand( P, 1 );
       
    [val ind_p ]  = max( p );
    
    ind_n = randi( N ); 
    
    if C( ind_p, ind_n ) ~= 0, continue, end 
    
    C( ind_p, ind_n ) = rand; 

    Sum_2(ind_p) = Sum_2(ind_p)+1;
    
    k=k+1; %fprintf( '%d\n', k ); 
end 
ind = find( sum( C, 1 ) == 0 );
for i=1:numel(ind)
    C( randi( P ), ind(i) ) = rand; 
end 