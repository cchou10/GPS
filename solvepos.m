% solvepos.m	(actual file name: solvepos.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this function calculates the position of the observation station
%
% input: 'ephem' 4-by-24 matrix, each of whose rows contains 
%                orbital ephemerides and other navigation message
%                data for a given satellite.

%					< see formatdata.m for description >
%
%		'pseudoR' 4-by-1 vector which contains pseudo-ranges (meters) for
%		          the satellites being used in the navigational solution
%						[ pr; pr; pr; pr]
%
%                 these must be pseudoranges from the same 4 satellites
%                 whose ephemerides and clock correction data are
%                 contained in the corresponding rows of ephem.
%
%		'guess'   3-by-1 vector containing an initial position guess in
%		          ECEF coordinates (meters).
%
%		'gpsTime' variable which contains the GPS time according to the
%                 receiver clock (seconds) for which a navigational  
%		          solution will be found.  This is the receiver
%                 clock time at which the 4 pseudoranges in pseudoR
%                 have been measured.
%
% output: 'posOBS'  1-by-5 row vector which contains the GPS receiver
%                   clock time (seconds), ECEF coordiates of the 
%                   navigation solution (meters), and the receiver clock  
%                   offset at that GPS time (seconds)
%						[ gpsTime, ECEFx, ECEFy, ECEFz, recCO]
%                   note that the true reception time corresponding
%                   to receiver clock time gpsTime is gpsTime_true = 
%                   gpsTime - recCO.
%
function posOBS = solvepos(ephem,pseudoR,guess,gpsTime)
% define physical constants
constant;
% compute the transmitter clock times associated with the receiver
% clock time in gpsTime and with the pseudoranges in pseudoR.
% also, ensure that the pseudoranges are stored in a column vector.
pseudoR_colvec = pseudoR(:);
t = gpsTime - pseudoR_colvec*(1/c);
% compute the ECEF positions of the 4 satellites at their transmission
% times along with their true transmission times and their transmitter
% clock errors.
satLocClock = findsatclock(ephem,t);



%find the transmitter clock error timer derivative and pass
%it as one of the output arguments
ccdot=satLocClock(:,7);



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
% initialize the guesses of the receiver position and of the
% receiver clock error.  also initialize a variable that holds
% the receiver clock error multiplied by the speed of light.
obsPos = guess(:);
recCO = min(deltatR) - 0.067;
c_recCO = c*recCO;
% initialize iteration counter.
iters=0;
% solve for position iteratively until solution is within an
% acceptable error
stop = 0;
while (stop == 0)
  iters=iters+1;
  % rotate satellite position vectors into ECEF reference frame
  % of the current guess of the true reception time.
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
  % compute the partial derivatives of the ranges with respect to the
  % elements of obsPos.
  oorange = range.^(-1);
  ax = - deltaXsatrcvr.*oorange;
  ay = - deltaYsatrcvr.*oorange;
  az = - deltaZsatrcvr.*oorange;
  % compute the partial derivatives of the ranges with respect recCO.
  arecCO = OmegaE*((-deltaXsatrcvr.*satY + deltaYsatrcvr.*satX).*oorange);
  %  form the vector 'l' and the matrix 'A'
  l = pseudoR_corrected - range + arecCO*recCO;
  A = [ax,ay,az,(ones(4,1) + arecCO*(1/c))];  
  % solve for "x", which contains deltaX, deltaY, deltaZ, and 
  % the new c*recCO
  x = A\l;
  % Update 'obsPos' by adding the first three elements of 'x' 
  % to the current guess
  obsPos = obsPos + x(1:3,1);
  % Update the receiver clock error estimate.  Also, compute
  % the change in this quantity multiplied by the speed of light.
  delta_c_recCO = x(4,1) - c_recCO;
  c_recCO = x(4,1);
  recCO = c_recCO/c;
  % check to see if the changes are small enough.  If they are,
  % then stop the iteration by setting the
  % 'stop' flag to 1; else, iterate again.
  normerror = norm([x(1:3,1);delta_c_recCO]);
  if ((normerror < 1.e-06) || (iters > 10)) %1.e-06 m change limit to stop
    stop = 1;
  end
  iters;
end
% create the output matrix 'posOBS' and return
posOBS = [gpsTime,obsPos',recCO, ccdot, range];
return;