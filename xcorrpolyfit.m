% xcorrpolyfit.m     (acutal file name: xcorrpolyfit.m)
%
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This function fits a polynomial to GPS pseudorange error data, 
% prerrhist, subtracts that polynomial from the data to create
% residual errors that have low-frequency trends removed,
% prresidhist, and compute the autocorrelation function of
% the resulting residuals time history.
%
%   inputs:   'prerrhist'       the N-by-1 vector pseudorange error 
%                               data from which to compute the   
%                               autocorrelation after subtracting out a 
%                               polynomial.  These should be sampled
%                               at roughly a constant frequency so that
%                               the time between samples prerrhist(k,1)
%                               and prerrhist(k+1,1) is independent of
%                               k or nearly independent of k.
%             'npolyfit'        The order of the polynomial that will be
%                               fit to the prerrhist data, fit as a  
%                               function of the row index of prerrhist.  
%                               Note, if npolyfit = 0, then a simple  
%                               computation of the mean and
%                               subtraction of it occurs to create
%                               the residuals whose autocorrelation 
%                               function is calculated.
%
%   outputs:  'Corrvec'         the Mlags-by-1 vector of normalized
%                               autocorelation values of prerrresidhist
%                               computed at the corresponding lags
%                               in mlagsvec.
%             'mlagsvec'        = (-round((N-1)/2):round((N-1)/2))', 
%                               the Mlags-by-1 vector of sample count 
%                               lags at which the autocorrelation 
%                               values in corrvec have been 
%                               calculated.  Note that Mlags = ...
%                               2*round((N-1)/2) + 1, that
%                               mlagsvec((round((N-1)/2) + 1),1) = 0, and 
%                               that corrvec((round((N-1)/2) + 1),1) = 1.
%             'polyprerr'       the 1-by-(npolyfit+1) vector of polynomial
%                               fit coefficients (meters) such that the 
%                               polynomial model of prerrhist is 
%                               prerrhistpoly = polyval(polyprerr,dnvec)
%                               where dnvec = (1:N)' - 0.5*(N+1)
%             'prerrresidhist'  = prerrhist - polyval(polyprerr,dnvec)
%                               where dnvec = (1:N)' - 0.5*(N+1), the
%                               N-by-1 vector of pseudorange
%                               residuals after subtraction of
%                               the polynomial in polyprerr.
%

function [ Corrvec, mlagsvec, polyprerr, prerrresidhist] = ...
                xcorrpolyfit(prerrhist,npolyfit)
%

% determine the number of points.
  prerrhistdum = prerrhist(:);
  N = size(prerrhistdum,1);
% do the polynomial fit and subtract it from the data to create the
% residuals.  Use a normalized version of dnvec in the polynomial
% fit in order to avoid numerical conditioning problems.
  dnvec = (1:N)' - 0.5*(N+1);
  dnvecmax = max(abs(dnvec));
  oodnvecmax = 1/dnvecmax;
  dnvec_nrmlzd = dnvec*oodnvecmax;
  polyprerr = polyfit(dnvec_nrmlzd,prerrhistdum,npolyfit);
  polyprerr = polyprerr.*((oodnvecmax*ones(1,(npolyfit+1))).^(npolyfit:-1:0));
  prerrhistpoly = polyval(polyprerr,dnvec);
  prerrresidhist = prerrhistdum - prerrhistpoly;
% determine the number of lags for the autocorrelation and set up
% the vector of lags.
  Mlagsm1o2 = round((N-1)/2);
  mlagsvec = (-Mlagsm1o2:Mlagsm1o2)';
% compute the unnormalized cross-correlation using the circular
% correlation function. xcorr is a special function that
% uses FFT techniques to rapidly compute the unnormalized sum in
% Eq. (8.7c) at many offset values.
  Corrvecun = xcorr(prerrresidhist,prerrresidhist,Mlagsm1o2,'none');
% normalize the correlations. note that D_reverse(k,1) = ...
% sqrt((prerrresidhist(1,1)^2 + ... + prerrresidhist(k,1)^2)*...
% (prerrresidhist(N-k+1,1)^2 + ... + prerrresidhist(N,1)^2)).  thus, it 
% contains the Dm coefficients of Eq. (8.7d) with Dm = D_reverse(N-m,1),
% hence the suffix "reverse" in its variable name.
  prerrresidhistr_reverse = prerrresidhist(N:-1:1,1);
  D_reverse = sqrt(cumsum(prerrresidhist.^2).*...
                  cumsum(prerrresidhistr_reverse.^2));
  Dvec = [D_reverse((N-Mlagsm1o2):N,1);D_reverse((N-1):-1:(N-Mlagsm1o2),1)];
  Corrvec = Corrvecun./Dvec;