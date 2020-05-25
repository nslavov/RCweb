function [B hat_B IND_OFF ind_corr] = SpAce_Path( V, P, Iter, p, Tresh, C_or_B ) 

% Sparse Ace (SpAce) for sparse factor analysis. Records the full path  
% Solves: hat_B = min_{B} ||VB||_0  so that ||B||_F = 1
%
%[B hat_B IND_OFF ind_corr] = SpAce_Path( V, P, Iter, p, Tresh, C_or_B )
%
% Input Parameters: Only V is required 
%
% V -- Input matrix. When SpAce is used with RCweb
%      for net inference it is the right singular matrix of the data  
% P -- Dimension of B; Number of hidden variables in RCweb
% Iter -- Maximum number of Iterations 
% p -- Number of hat_B columns to be infered. Default p = P; 
% Thresh -- Termination threshold for the eigenvalue of inv(V(on,:)*'V(on,:))
% C_or_B -- Used only for evaluation of net inference from simulated data.
%           Gold standard B or C (adjecancy matrix) for the simulated net 
%
%
% Output Parameters:
% hat_B -- optimal B
% hat_B -- Cell arrays containing smlp number of samples from b_i  
% IND_OFF -- Cell containg the set of indecies that are 
%             non-zero in each column of ||V hat_B||_0  
%ind_corr -- Used only for evaluation of net inference from simulated data.
%            Correlations (and associated permutation indecies) between   
%            C_or_B and hat_B or the infered adjecany matrix (C=V*hat_B)
%
% Nikolai Slavov, April.2009


% == Parameter Definitions and setting Defaults ==
if nargin < 2 || isempty(P), P = size(V,2);      end
if nargin < 3 || isempty(Iter), Iter = 1e2;      end
if nargin < 4 || isempty(p), p = P;              end 
if nargin < 5 || isempty(Tresh), Tresh = 1e3;    end 
if nargin >=6
    if  size(C_or_B,2) == P,         TEST = 1;   end  % B
    if  size(C_or_B,2) == size(V,1), TEST = 2;   end  % C
        ind_corr = zeros( p, 2 );
else TEST = [];
end 
         
defl = round(max( 20, round(P/25) )); % Number of Iter with the deflated matrix 

% == variable Initiation ==
hat_B = zeros( P, p );
B = cell( p, 1 );
IND_OFF = cell( p, 1 );

smlp = 100; %Number of samples 
SMPL = round( linspace( 1, Iter,  smlp ) );  
b = zeros( P, smlp ); 
Sing = zeros( Iter, p ); % Array to store the singular values 


% == Initialization for a Matrix W=U whose columns are orthogonal 
[val ind_off] = max( sum(abs(V),2) ); %% Chose an index to remove
IND_off(1) = ind_off; 

iK = eye( P ); 
ik = inv_rank1_red( iK, V(ind_off,:)' );  v = rand(P,1); 

% == Main Loop
for I = 1:p   
    %fprintf( 'Iteration: %d\t%1.0f s\n', I, (toc-t) );  t=toc;
     
    smpl = 1;
    if I>1, %//== BLOCK 1 ==//        
        IND_off = [];
        [u lam] = pm( Defl'*Defl, [], 1e3 );  %lam=2*lam;    
        iKu = iK*u;
        sk = 1/( 1/lam + u'*iK*u );
        iK = iK - sk*iKu*iKu';
        ik = iK; 
    end
    
    for i=1:Iter % CYCLE :: Breaks Inside when the smallest singular value --> 0   

        [v lam] = pm( ik, v, 5e2 );   %fprintf( 'Point 1: %1.2f s\n', toc );        
          
        Sing(i,I) = lam; 
        if i==SMPL(smpl)
            b(:,smpl) = v;  smpl=smpl+1;
        end 
        
        if abs(lam)>Tresh,  break, end % Breaks if W_{on} --> singular 1e+7 
        

        u = V*v;
        u(IND_off) = 0;
        [val ind_off] = max( abs(u) );  
        IND_off = [IND_off; ind_off]; %#ok<*AGROW>
        ik = inv_rank1_red( ik, V(ind_off,:)' );
        
        
        if I>1 && i==defl   %//== BLOCK 2 ==//
             rows_off = V( IND_off, : );
             ik = inv_up( eye( P ),  rows_off',  eye(numel(IND_off)),...
                                     rows_off,    -1     );
        end 
    end 
    Defl = (V*v) * v'; %W1 = W1 + Defl;

    
    B{I,1} = b(:,1:smpl-1);
    hat_B(:,I) = v;
    IND_OFF{I,1} = IND_off;
    
     
    switch TEST
        case { 1, 'B' } 
            cor = corr( C_or_B, hat_B(:,I) ); 
            [ind_corr(I,1) ind_corr(I,2)] = max( abs(cor) ); 
            hat_B(:,I) = hat_B(:,I)* norm(C_or_B(:,ind))*sign(cor(ind));        
        case { 2, 'C' }
            cor = corr( C_or_B', V*v );
            [ind_corr(I,1) ind_corr(I,2)] = max( abs(cor) );
    end 
    if ~isempty(TEST)
           fprintf( '%d) Max Corr: %1.2f : %d\n', I, ind_corr(I,:) ); 
    end    
end
