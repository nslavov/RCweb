function [u lam_2 i] = pm( x, u, Iter, Thresh )


if nargin<2  || isempty(u), u = rand( size(x,2), 1 ); end  
if nargin<3  || isempty(Iter), Iter = 3e2;            end 
if nargin<4, Thresh=1e-7;               end 

u_norm = norm( u );
lam_1 = -Inf;  

for i=1:Iter
    
   %old = u;
    old_norm = u_norm;
     
    u = x*u;
    
    u_norm = norm( u ); 
    lam_2 = u_norm / old_norm;
    
    pr = abs(lam_2-lam_1)/lam_2;
    
    if i>1 &&  pr<Thresh,  break; end 
    
    lam_1 = lam_2;
    
    if u_norm >1e30, u = (1e-10/u_norm)*u;  end 
end

u = u * (1/u_norm);

if isnan(lam_2)
    fprintf( 'The pm.m returned NaN !\n' );  
end 

if i==Iter && pr>1e2*Thresh
   fprintf( '%d iterations were not enough to converge: %1.0e !\n', i, pr ); 
end  
