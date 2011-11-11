function [ Ke_, fue_ ] = cor_compute_stiffness_warp(Re, Ke, fue, params)
% Copyright 2011, Kenny Erleben

Ke_  = Ke;
fue_ = fue;

if( params.warp )
  
  E = length(Re(1,:)) / 3;
  
  for e=1:E
    Re_  = Re(:,( 3*(e-1)+1):( 3*e));
    R    = blkdiag(Re_,Re_,Re_,Re_);
    
    Ke_(:,(12*(e-1)+1):(12*e))  = R * Ke(:,(12*(e-1)+1):(12*e)) * R';
    fue_(:,e)                   = R*fue(:,e);
    
  end
  
end

end