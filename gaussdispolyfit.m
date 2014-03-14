% gaussdispolyfit.m     (acutal file name: gaussdispolyfit.m)
%
%  Copyright (c) 2000 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This function fits a polynomial to GPS pseudorange error data, 
% prerrhist, subtracts that polynomial from the data to create
% residual errors that have low-frequency trends removed,
% prresidhist, computes the standard deviaton of the resulting 
% residuals, sigmapr, creates a histogram of the residuals by  
% counting the data in bins 0.5*sigmapr wide, from -3.0*sigmapr to 
% 3.0*sigmapr around mean(prresidhist) -- this mean should be zero,
% computing the expected numbers in the bins based on the Gaussian
% assumption, plotting this raw histogram and expected number
% of counts data, and computing the Pearson chi-squared chi-squared
% test statistic.
%
%   inputs:   'prerrhist'     the N-by-1 vector pseudorange error 
%                             data from which to create   
%                             the histogram after subtracting out a 
%                             polynomial.
%             'npolyfit'      The order of the polynomial that will be
%                             fit to the prerrhist data, fit as a  
%                             function of the row index of prerrhist.  
%                             Note, if npolyfit = 0, then a simple  
%                             computation of the mean and
%                             subtraction of it occurs to create
%                             the residuals that are binned into
%                             a histogram.
%
%   outputs:  'sigmapr'       standard deviation of the residual pseudo
%                             ranges after subtracking of the polynomial
%                             (meters).
%             'polyprerr'     the 1-by-(npolyfit+1) vector of polynomial
%                             fit coefficients (meters) such that the 
%                             polynomial model of prerrhist is 
%                             prerrhistpoly = polyval(polyprerr,dnvec)
%                             where dnvec = (1:N)' - 0.5*(N+1)
%             'Ovec'          the 12-by-1 histogram vector that contains
%                             the counts of the histogram in the
%                             bins whose lower and upper bounds,
%                             expressed as standard deviations from
%                             the mean, are given in binlohihisto.
%                             Note: in order to lie in bin k, i.e.,
%                             in order to be counted in the number
%                             Ovec(k,1), it must be greather than
%                             or equal to the lower bin boundary
%                             in binlohihisto(k,1) and it must be
%                             less than the upper bin boundary
%                             in binlohihisto(k,1).
%             'Evec'          the 12-by-1 vector that is like Ovec,
%                             except that it is the expected number
%                             of counts given the Gaussian assumption
%                             for the residuals in prresidhist with
%                             zero mean and a standard deviaton of
%                             sigmapr.
%             'M'             = 12, the number of histogram bins.
%             'bin_pos'       the 12-by-1 vector of mid-points of
%                             prresidhist bins = (-2.75 + (0:11)'*0.5)*...
%                             sigmapr (meters).
%             'N'             the total number of data points.
%             'chi2Pearson'   The Pearson chi-squared test statistic
%                             for the goodness of the Gaussian model
%                             (non-dimensional)
%             'ndegreechi2'   = M - 3 - npolyfit, the modeled
%                             degree of the Pearson chi-squared test
%                             statistic.
%             'chi2Pearson95' = chi2inv(.95,ndegreechi2), the 95%
%                             threshold that chi2Pearson should
%                             be below in 95% of the cases where the
%                             distribution really is Gaussian.
%             'chi2Pearson99' = chi2inv(.99,ndegreechi2), the 99%
%                             threshold that chi2Pearson should
%                             be below in 99% of the cases where the
%                             distribution really is Gaussian.
%
%             Note: chi2Pearson95 and chi2Pearson99 will be empty 
%             on output of the function chi2inv is not available
%             due to there not being a statistics toolbox.
%

function [ sigmapr, polyprerr, Ovec, Evec, M, bin_pos, ...
                   N, chi2Pearson, ndegreechi2, chi2Pearson95, ...
                   chi2Pearson99] = ...
                gaussdispolyfit(prerrhist,npolyfit)
%

% use hist3sigpolyfit.m to compute the histogram.
  [ mprresid, sigmapr, polyprerr, Ovec, M, binlohihisto, ...
                   bin_pos, N ] = ...
                hist3sigpolyfit(prerrhist,npolyfit);
% compute the theoretical values for the histogram
  Prgauss = [0.004860;...
             0.016540;...
             0.044057;...
             0.091848;...
             0.149882;...
             0.191462;...
             0.191462;...
             0.149882;...
             0.091848;...
             0.044057;...
             0.016540;...
             0.004860];
  Evec = N*Prgauss;
% compute the Pearson chi-squared test statistic and its theoretical
% degree
  chi2Pearson = sum(((Ovec - Evec).^2)./Evec);
  ndegreechi2 = M - 3 - npolyfit;
% compute the chi-squared threshold values if the function chi2inv.m 
% exists on a path (exist('chi2inv') == 2) or if it is a built-in
% function (exist('chi2inv') == 5). otherwise, do not calculate
% thresholds, but display a warning.
  exist_chi2inv = exist('chi2inv');
  if (exist_chi2inv == 2) | (exist_chi2inv == 5)
     chi2Pearson95 = chi2inv(0.95,ndegreechi2);
     chi2Pearson99 = chi2inv(0.99,ndegreechi2);
  else
     chi2Pearson95 = [];
     chi2Pearson99 = [];
     disp('Warning in gaussdispolyfit.m: The chi-squared threshold')
     disp(' outputs chi2Pearson95 and chi2Pearson99 could not be')
     disp(' computed because the inverse chi-squared probability')
     disp(' function chi2inv.m is not available.  The statistics')
     disp(' toolbox is probably not installed on this computer.')
  end
% plot the histogram and its theoretical Gaussian values.
  figure
  hold off
  bar(bin_pos,Ovec);
  grid
  hold on
  set(get(gcf,'CurrentAxes'),'FontSize',16);
  plot(bin_pos,Evec,'r*');
  hold off
  xlabel('pseudorange residual')
  ylabel('histogram bin count')
  legend('Experimental','Theoretical Gaussian')
  title(['Histogram of Pseudo-range Residuals after ',int2str(npolyfit),...
         '-order polynomial detrending'])
  hold off
  return