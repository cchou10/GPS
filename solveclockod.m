% solveclockod.m	(actual file name: solveclockod.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this function calculates the receiver clock error for 
% a receiver at a known location.  It also uses the known
% location, the solved-for receiver clock error, and the
% pseudoranges to compute pseudorange errors for
% all signals.
%
% This function has the option to input an alternate
% set of pseudoranges for which errors are also computed even
% though these alternate pseudoranges are not used for
% the clock offset solution.  The reason for allowing
% a second set of pseudoranges for the same satellites
% is that the pseudoranges used for the actual clock error
% solution might be beat-carrier-phase-based pseudoranges
% that equal beat carrier phases multiplied by nominal
% carrier wave length with biases removed to have the same
% average over a data span as the raw pseudoranges.  The 
% second set of pseudoranges can be the raw pseudoranges.
% The corresponding pseudorange errors will allow the study
% of both beat-carrier-phase measurement error statistics
% and pseudorange measurement error statistics.
%
% This function also has an input vector that allows for
% the possibility that only a subset of the input pseudoranges
% will be used for determining the clock offset solution.
% This feature is useful if not all available signals are
% available for the entire time span of interest.  It is best
% to use the same set of satellites for determining the
% clock error over the entire time span of interest.  This
% feature allows for that possibility.
%
% This function also computes the azimuths and elevations of the 
% satellites.
%
%
% input: 'ephem'    Nsats-by-24 matrix, each of whose rows contains 
%                   orbital ephemerides and other navigation message
%                   data for a given satellite.
%
%					  < see formatdata.m for description >
%
%		 'pseudoR'  Nsats-by-1 vector which contains pseudo-ranges
%                   (meters) for the satellites being used in the 
%		            clock error solution and for use in the
%                   pseudorange error calculations.
%
%						[ pr; pr; pr; ...; pr]
%
%                   These must be pseudoranges from the same satellites
%                   whose ephemerides and clock correction data are
%                   contained in the corresponding rows of ephem.
%                   Note that, depending on the values in 'iflagpR',
%                   not all the entries of pseudoR may be used
%                   for clock error calculations.
%
%        'iflagpR'  Nsats-by-1 vector of flags that indicate which
%                   of the pseudoranges in pseudoR are to be used
%                   to compute the receiver clock error.  The
%                   used set of pseudoranges is that for which
%                   the corresponding entries of 'iflagpR' equal 1.
%
%		 'obsPos'   3-by-1 vector containing the known position of
%                   the receiver in ECEF coordinates (meters).
%
%		 'gpsTime'  variable which contains the GPS time according to the
%                   receiver clock (seconds) for which the receiver  
%		            clock error will be found.  This is the receiver
%                   clock time at which the Nsats pseudoranges in pseudoR
%                   have been measured.

%		 'pseudoR2' Nsats-by-1 vector which contains alternate
%                   pseudo-ranges (meters) for the satellites 
%                   in order to calculate alternate pseudorange
%                   errors.
%
%						[ pr2; pr2; pr2; ...; pr2]
%
%                   These must be pseudoranges from the same satellites
%                   whose ephemerides and clock correction data are
%                   contained in the corresponding rows of ephem.
%                   These pseudoranges are not used to compute the
%                   receiver clock error.  Rather, they
%                   are used only to compute alternate pseudorange
%                   errors.  This is an optional input.  If it is
%                   not input, then the alternate pseudorange error
%                   output pseudoerror2 will be an empty array.
%                   Typically the pseudoranges in pseudoR are the
%                   ones with the lowest measurement errors, e.g.,
%                   carrier-smoothed pseudoranges, while the ones 
%                   in pseudoR2 are ones with larger measurement
%                   errors, e.g., raw pseudoranges.
%
% outputs: 'recCO'          The scalar receiver clock offset (seconds)
%					        Note that the true reception time
%                           corresponding to receiver clock time 
%                           gpsTime is gpsTime_true = gpsTime - recCO.
%
%          'pseudoerror'    The Nsats-by-1 vector of pseudorange errors 
%                           (meters) associated with the pseudoranges
%                           in pseudoR, with the known receiver position
%                           in obsPos, and with the receiver clock error
%                           estimate in recCO.
%
%						      [ prerr; prerr; prerr; ...; prerr]
%
%          'pseudoerror2'   The Nsats-by-1 vector of pseudorange errors 
%                           (meters) associated with the pseudoranges
%                           in pseudoR2, the with known receiver position
%                           in obsPos, and with the receiver clock error
%                           estimate in recCO.
%
%						      [ prerr2; prerr2; prerr2; ...; prerr2]
%
%                           Note: this is an empty array on output
%                           if pseudoR2 is empty or omitted on
%                           input.
%
%          'azimvec'        The Nsats-by-1 vector of the azimuths
%                           (deg from North, with East equal to 90 deg)
%                           associated with the satellites who SVIDs
%                           are contained in ephem(:,1).
%
%          'elevvec'        The Nsats-by-1 vector of the elevations
%                           (deg) associated with the satellites who 
%                           SVIDs are contained in ephem(:,1).
%
function [recCO,pseudoerror,pseudoerror2,azimvec,elevvec] = ...
           solveclockod(ephem,pseudoR,iflagpR,obsPos,gpsTime,pseudoR2)
% define physical constants
constant;
% get the number of satellites to use in the pseudorange calculations.
Nsats = size(ephem,1);
% get the indices of and the number of satellites to use in the 
% receiver clock error solution. stop with an error if there are
% too few satellites.
iuseforrecCO = find(iflagpR == 1);
Nsats_sht = size(iuseforrecCO,1);
if Nsats_sht < 1
   error(['Error in solveclockod.m: fewer than 1 satellite',...
          ' pseudorange entered for use in clock error',...
          ' calulation, i.e., no 1 flag values were',...
          ' entered in iflagpR.'])
end
% compute the transmitter clock times associated with the receiver
% clock time in gpsTime and with the pseudoranges in pseudoR.
% also, ensure that the pseudoranges are stored in a column vector.
pseudoR_colvec = pseudoR(:);
t = gpsTime - pseudoR_colvec*(1/c);
% compute the ECEF positions of the Nsats satellites at their transmission
% times along with their true transmission times and their transmitter
% clock errors.
satLocClock = findsatclock(ephem,t);
% compute the time difference vector deltatR such that
% deltprop = deltatR - recCO is the vector of signal propagation
% delays, in seconds, if recCO is the receiver clock error.
deltatR = gpsTime - satLocClock(:,2);
% compute pseudoranges that have been corrected for the transmitter
% clock errors.
pseudoCorr = satLocClock(:,6)*c;
pseudoR_corrected = pseudoR_colvec + pseudoCorr;
% set up the GPS satellite ECEF X and Y positions in the
% transmission time coordinate frames.
satX_Trans = satLocClock(:,3);
satY_Trans = satLocClock(:,4);
% set up the GPS satellite ECEF Z positions.
satZ = satLocClock(:,5);
% set up shortened nominal transmission delay, pseudorange
% and position vectors for use in the clock error solution.
deltatR_sht = deltatR(iuseforrecCO,1);
pseudoR_corrected_sht = pseudoR_corrected(iuseforrecCO,1);
satX_Trans_sht = satX_Trans(iuseforrecCO,1);
satY_Trans_sht = satY_Trans(iuseforrecCO,1);
satZ_sht = satZ(iuseforrecCO,1);
% initialize the guess of the
% receiver clock error.  also initialize a variable that holds
% the receiver clock error multiplied by the speed of light.
recCO = min(deltatR_sht) - 0.067;
c_recCO = c*recCO;
% initialize iteration counter.
iters=0;
% solve for clock error iteratively until solution is within an
% acceptable error
stop = 0;
while (stop == 0)
  iters=iters+1;
  % rotate satellite position vectors into ECEF reference frame
  % of the current guess of the true reception time.
  deltheta_sht = OmegaE*(deltatR_sht - recCO);
  cos_deltheta_sht = cos(deltheta_sht);
  sin_deltheta_sht = sin(deltheta_sht);
  satX_sht = cos_deltheta_sht.*satX_Trans_sht + ...
               sin_deltheta_sht.*satY_Trans_sht;
  satY_sht = -sin_deltheta_sht.*satX_Trans_sht + ...
               cos_deltheta_sht.*satY_Trans_sht;
  % compute the ranges to the satellites.
  deltaXsatrcvr_sht = satX_sht - obsPos(1,1);
  deltaYsatrcvr_sht = satY_sht - obsPos(2,1);
  deltaZsatrcvr_sht = satZ_sht - obsPos(3,1);
  range_sht = sqrt(deltaXsatrcvr_sht.^2 + deltaYsatrcvr_sht.^2 + ...
                   deltaZsatrcvr_sht.^2);
  % compute the inverse range.
  oorange_sht = range_sht.^(-1);
  % compute the partial derivatives of the ranges with respect recCO.
  arecCO_sht = OmegaE*((-deltaXsatrcvr_sht.*satY_sht + ...
                        deltaYsatrcvr_sht.*satX_sht).*oorange_sht);
  % form the shortened versions of the vector 'l' and the 
  % "matrix" 'A1col'.  The latter "matrix" includes only the final
  % column of the usual 'A' matrix because this function assumes
  % the positioin in obsPos is correct, and it solves only for the
  % receiver clock error, which is the unknown that multiplies
  % the last column of the original 'A' matrix.
  l_sht = pseudoR_corrected_sht - range_sht + arecCO_sht*recCO;
  A1col_sht = (ones(Nsats_sht,1) + arecCO_sht*(1/c));  
  % solve the possibly over-determined problem for the new c*recCO 
  x4 = A1col_sht\l_sht;
  % Update the receiver clock error estimate.  Also, compute
  % the change in this quantity multiplied by the speed of light.
  delta_c_recCO = x4 - c_recCO;
  c_recCO = x4;
  recCO = c_recCO/c;
  % check to see if the change is small enough.  If it is,
  % then stop the iteration by setting the
  % 'stop' flag to 1; else, iterate again.
  normerror = abs(delta_c_recCO);
  if ((normerror < 1.e-06) || (iters > 10)) %1.e-06 m change limit to stop
    stop = 1;
  end
  iters;
end
% do the full pseudorange error computations for all pseudoranges.
deltheta = OmegaE*(deltatR - recCO);
cos_deltheta = cos(deltheta);
sin_deltheta = sin(deltheta);
satX = cos_deltheta.*satX_Trans + sin_deltheta.*satY_Trans;
satY = -sin_deltheta.*satX_Trans + cos_deltheta.*satY_Trans;
% compute the ranges to the satellites.
deltaXsatrcvr = satX - obsPos(1,1);
deltaYsatrcvr = satY - obsPos(2,1);
deltaZsatrcvr = satZ - obsPos(3,1);
range = sqrt(deltaXsatrcvr.^2 + deltaYsatrcvr.^2 + deltaZsatrcvr.^2);
% form psdudorange error vector
rangepc_recCO = range + c_recCO;
pseudoerror = pseudoR_corrected - rangepc_recCO;
% form the alternate pseudorange error vector if alternate
% pseudoranges have been entered in pseudoR2.  otherwise,
% output an empty array for pseudoerror2.
pseudoerror2 = [];
if nargin >= 6
   if size(pseudoR2,1) == Nsats
      pseudoerror2 = pseudoR2 + pseudoCorr - rangepc_recCO;
   end
end
% compute the satellite azimuths and elevations and place them in
% the proper output arrays.
el_az = elevazim([ephem(:,1),t,satX,satY,satZ],obsPos');
azimvec = el_az(:,4);
elevvec = el_az(:,3);
% send a warning message to the display if any elevations are
% below 10 deg among the set of satellites used to
% determine the receiver clock error.
idum = find(elevvec(iuseforrecCO,1) < 10);
ndum = size(idum,1);
if ndum >= 1
   idum  = iuseforrecCO(idum,1);
   for jjdum = 1:ndum
      jj = idum(jjdum,1);
      SVjj = ephem(jj,1);
      if SVjj < 10
          SVjjtext = ['0',int2str(SVjj)];
      else
          SVjjtext = int2str(SVjj);
      end
      disp(['Warning in solveclockod.m: PRN',SVjjtext,...
            ' has elevation ,',num2str(elevvec(jj,1)),...
            ' deg, i.e., below the recommended 10 deg mask.'])
   end
end
return;