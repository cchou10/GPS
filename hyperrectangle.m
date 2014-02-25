% hyperrectangle.m	(actual file name: hyperrectangle.m)
%
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
%	< USE THE REAL-VALUED ESTIMATE, ITS SIGMAS, AND A NON-DIMENSIONAL 
%     SIGMA FACTOR IN ORDER TO DEVELOP A LIST OF ALL INTEGER-
%     VALUED VECTOR LATTICE POINTS IN THE HYPER-RECTANGLE 
%     CENTERED AT THE REAL-VALUED ESTIMATE. >
%
% this CDGPS utility function computes all of the double-differenced
% beat carrier phase ambiguity integer points that lie inside or on
% a given hyper-rectangle.  This hyper-rectangle is centered at the
% real-valued estimate of the ambiguity vector, and the halfwidths
% of its sides are defined by a prescribed number of theoretical 
% standard deviations of its real-valued estimate's components.
%
% input:  'graddeltabeta_r'        The (N-1)-dimensional vector of
%                                  real-valued estimates of the
%                                  doubled-differenced beat carrier 
%                                  phase ambiguities.
%
%         'sigma_graddeltabeta_r'  The (N-1)-dimensional vector of
%                                  theoretical standard deviations
%                                  of the real-valued estimates of
%                                  the doubled-differenced beat  
%                                  carrier phase ambiguities.
%
%         'aeta'                   The positive scalar number of
%                                  standard deviations in
%                                  sigma_graddeltabeta_r that will
%                                  be used to define the half
%                                  width of the hyper-rectangle
%                                  in each direction.  The hyper-
%                                  rectangle will be centered on 
%                                  the real-valued estimate in 
%                                  graddeltabeta_r.
%
%                                  Note: The MAE 4150/ECE 4150
%                                  formulas of sigma_graddeltabeta_r
%                                  and aeta depend on knowledge of
%                                  the un-differenced carrier-phase
%                                  measurement error standard 
%                                  deviation sigmaphi.  In practice,
%                                  any arbitrary positive value for
%                                  sigmaphi can be assumed in computing
%                                  both sigma_graddeltabeta_r and
%                                  aeta without changing the results
%                                  of the present function.  The
%                                  most convenient arbitrary value
%                                  is sigmaphi = 1.
%
% output: graddeltabeta_inrect'    The (N-1)-by-Ninrectangle matrix,
%                                  each of whose columns contains one
%                                  of the vector integer-ambiguity
%                                  lattice points in the hyper-rectangle
%                                  that is defined by the present
%                                  function's inputs.  This array
%                                  contains all of the integer-valued
%                                  vector lattice points in the 
%                                  defined hyper-rectangle.
%
%	< USE THE REAL-VALUED ESTIMATE, ITS SIGMAS, AND A NON-DIMENSIONAL 
%     SIGMA FACTOR IN ORDER TO DEVELOP A LIST OF ALL INTEGER-
%     VALUED VECTOR LATTICE POINTS IN THE HYPER-RECTANGLE 
%     CENTERED AT THE REAL-VALUED ESTIMATE. >
%
function graddeltabeta_inrect = hyperrectangle(graddeltabeta_r,...
                                               sigma_graddeltabeta_r,...
                                               aeta)
% determine the lower and upper limits of the range of integer
% ambiguities that will be searched.  make double-sure that each 
% range includes the corresponding rounded solution.
graddeltabeta_round = round(graddeltabeta_r);
graddeltabeta_lo = ceil(graddeltabeta_r - aeta*sigma_graddeltabeta_r);
graddeltabeta_lo = min([graddeltabeta_lo,graddeltabeta_round],[],2);
graddeltabeta_hi = floor(graddeltabeta_r + aeta*sigma_graddeltabeta_r);
graddeltabeta_hi = max([graddeltabeta_hi,graddeltabeta_round],[],2);
% compute all the possible integer vectors, one for each column of
% the array graddeltabeta_inrect
graddeltabeta_counts = graddeltabeta_hi - (graddeltabeta_lo - 1);
Ninrectangle = prod(graddeltabeta_counts);
Nm1 = size(graddeltabeta_r,1);
graddeltabeta_inrect = zeros(Nm1,Ninrectangle);
mrepeatjj = 1;
for jj = 1:Nm1
   idumvecjj = (graddeltabeta_lo(jj,1):graddeltabeta_hi(jj,1));
   ndumjj = graddeltabeta_counts(jj,1);
   idumvecjj = ones(mrepeatjj,1)*idumvecjj;
   idumvecjj = idumvecjj(:);
   mrepeatjj = mrepeatjj*ndumjj;
   mremainingjj = Ninrectangle/mrepeatjj;
   idumvecjj = idumvecjj*ones(1,mremainingjj);
   idumvecjj = idumvecjj(:);
   graddeltabeta_inrect(jj,:) = idumvecjj';
end
return;