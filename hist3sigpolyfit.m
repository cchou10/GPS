% hist3sigpolyfit.m	(acutal file name: hist3sigpolyfit.m)
%
%  Copyright (c) 2000 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This function fits a polynomial to random data, subtracts
% that polynomial from the data, computes the standard deviaton
% of the resulting residuals, and creates a histogram of the
% residuals by counting the data in bins
% 0.5*std(xresid) wide, from -3.0*std(xresid) to 3.0*std(xresid)
% around mean(xresid).
%
%   inputs:   'x'             the N-by-1 vector data from which to create   
%                             the histogram after subtracting out a 
%                             polynomial.
%             'npolyfit'      The order of the polynomial that will be
%                             fit to the x data, fit as a function of
%                             the row index of x.  Note, if npolyfit = 0,
%                             then a simple computation of the mean
%                             and subtraction of it occurs to create
%                             the residuals that are binned into
%                             a histogram.
%
%   outputs:  'mxresid'       mean of xresid = x - polyval(polyx,dnvec) 
%                             where dnvec = (1:N)' - 0.5*(N+1).
%                             This should be zero.
%			  'sxresid'       standard deviation of xresid.
%             'polyx'         the 1-by-(npolyfit+1) vector of polynomial
%                             fit coefficients such that the polynomial
%                             model of x is xpoly = polyval(polyx,dnvec)
%                             where dnvec = (1:N)' - 0.5*(N+1)
%             'histo'         the 12-by-1 histogram vector that contains
%                             the counts of the histogram in the
%                             bins whose lower and upper bounds,
%                             expressed as standard deviations from
%                             the mean, are given in binlohihisto.
%                             Note: in order to lie in  bin k, i.e.,
%                             in order to be counted in the number
%                             histo(k,1), it must be greather than
%                             or equal to the lower bin boundary
%                             in binlohihisto(k,1) and it must be
%                             less than the upper bin boundary
%                             in binlohihisto(k,1).
%             'n_bins'        = 12, the number of histogram bins.
%             'binlohihisto'  The 12-by-2 array of histogram bin boundaries
%                             of xresid.  It equals [(-3 + (0:11)'*0.5),...
%                             (-3 + (1:12)'*0.5)]*sxresid + mxresid.
%             'bin_pos'       the 12-by-1 vector of mid-points
%                             xresid bins = (-2.75 + (0:11)'*0.5)*...
%                             sxresid + mxresid.
%             'tot_pts'       = N, the total number of data points.
%             'xresid'        = x - polyval(polyx,dnvec) where 
%                             dnvec = (1:N)' - 0.5*(N+1), the residuals
%                             history.
%

function [ mxresid, sxresid, polyx, histo, n_bins, binlohihisto, ...
                   bin_pos, tot_pts, xresid] = hist3sigpolyfit(x,npolyfit)
%

% determine the number of points.
  xdum = x(:);
  N = size(xdum,1);
  tot_pts = N;
% do the polynomial fit and subtract it from the data to create the
% residuals.  Use a normalized version of dnvec in the polynomial
% fit in order to avoid numerical conditioning problems.
  dnvec = (1:N)' - 0.5*(N+1);
  dnvecmax = max(abs(dnvec));
  oodnvecmax = 1/dnvecmax;
  dnvec_nrmlzd = dnvec*oodnvecmax;
  polyx = polyfit(dnvec_nrmlzd,xdum,npolyfit);
  polyx = polyx.*((oodnvecmax*ones(1,(npolyfit+1))).^(npolyfit:-1:0));
  xpoly = polyval(polyx,dnvec);
  xresid = xdum - xpoly;
% compute the mean and standard deviation of the residuals.
  mxresid = mean(xresid);
  sxresid = sqrt(sum(xresid.^2)*(1/(N - (npolyfit + 1))));
% construct the bin boundaries and centers.
   binlohihisto = [(-3 + ((0:11)')*0.5),(-3 + ((1:12)')*0.5)]*sxresid + ...
                  mxresid;
   bin_pos = (-2.75 + ((0:11)')*0.5)*sxresid + mxresid;
   n_bins = size(bin_pos,1);
% count the number of points in each bin.
   histo = zeros(n_bins,1);
   for jj = 1:n_bins
      histo(jj,1) = sum((binlohihisto(jj,1) <= xresid) & ...
                        (xresid < binlohihisto(jj,2)));
   end
   return