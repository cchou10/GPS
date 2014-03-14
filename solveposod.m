% solveposod.m	(actual file name: solveposod.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% this function calculates the position of the observation station
% using 4 or more satellites.  It also uses the ionosphere
% and neutral atmosphere corrections.
%
% input: 'ephem'   Nsatsall-by-24 matrix, each of whose rows contains 
%                  orbital ephemerides and other navigation message
%                  data for a given satellite.
%
%					< see formatdata.m for description >
%
%		 'pseudoR'  Nsatsall-by-1 vector which contains pseudo-ranges
%                   (meters) for the satellites being used in the 
%		            navigational solution
%
%						[ pr; pr; pr; ...; pr]
%
%                   these must be pseudoranges from the same satellites
%                   whose ephemerides and clock correction data are
%                   contained in the corresponding rows of ephem.
%
%		 'guess'    3-by-1 vector containing an initial position guess in
%		            ECEF coordinates (meters).
%
%		 'gpsTime'  variable which contains the GPS time according to the
%                   receiver clock (seconds) for which a navigational  
%		            solution will be found.  This is the receiver clock 
%                   time at which the Nsatsall pseudoranges in pseudoR
%                   have been measured.
%
%        'ionParam' 4-by-2 matrix of the eight coefficients 
%                   used to correct for ionospheric delay
%                   (seconds in first row, seconds/semi-circle
%                   in second row, seconds/(semi-circle)^2 in
%                   third row, and seconds/(semi-circle)^3 in
%                   4th row.)
%
%                          [ alpha0 beta0
%                            alpha1 beta1
%                            alpha2 beta2
%                            alpha3 beta3 ]
%
%        'iflagion' Flag that tells whether to use the ionospheric
%                   delay corrections in the calculation of
%                   the navigation solution (iflagion = 1) or
%                   whether to neglect them (iflagion = 0).
%
%        'elevmask' The elevation mask angle (degrees) below
%                   which a given satellite's data will not be
%                   used.  A recommended value is something in
%                   the range 5-10 deg.
%
%		 'p'        Atmospheric pressure at receiver (millibar or,
%                   identically, hPa)
%
%        'TdegK'    Temperature at receiver (degrees Kelvin)
%
%        'hrel'     Relative humidity at receiver (fraction in 
%                   range 0 to 1)
%
%                   Note: the last three entries are optional.
%                   If they are not included, then the following 
%                   "standard" sea-level model will be used:
%
%                      p     = 1013.25 millibars
%                      TdegK =  288.15 deg K
%                      hrel  =    0.50
%
%                   These values of p and TdegK are from the
%                   standard atmosphere.  The hrel value is
%                   a made-up standard.  The actual data may be
%                   retrieved from a web site such as 
%                   http://www.wunderground.com using a meteorological
%                   station in the same town as the receiver.  This  
%                   web site includes history data for past days.   
%                   Note that there are 33.86 millbars per inch of 
%                   Mercury (in Hg), that TdegK = TdegC + 273.15,
%                   and that TdegC = (TdegF - 32)*(5/9).
%
%        'iflagna'  Flag that tells whether to use the neutral
%                   atmosphere delay corrections in the calculation of
%                   the navigation solution (iflagna = 1) or
%                   whether to neglect them (iflagna = 0).
%
% output: 'posOBS'  1-by-5 row vector which contains the GPS receiver
%                   clock time (seconds), ECEF coordiates of the 
%                   navigation solution (meters), and the receiver clock  
%                   offset at that GPS time (seconds)
%
%						[ gpsTime, ECEFx, ECEFy, ECEFz, recCO]
%
%                   note that the true reception time corresponding
%                   to receiver clock time gpsTime is gpsTime_true = 
%                   gpsTime - recCO.
%
%         'DOP'     5-by-1 column vector which contains the various 
%                   measures of dilution of precision
%
%                     DOP(1,1) geometrical dilution of precision
%                     DOP(2,1) positional dilution of precision
%                     DOP(3,1) time dilution of precision
%                     DOP(4,1) horizontal dilution of precision
%                     DOP(5,1) vertical dilution of precision
%
%         'el_az'   Nsatsall-by-4 matrix each of whose rows contain an
%                   SV id number, a GPS time (seconds), and the 
%                   elevation and azimuth look angles (degrees) to the 
%                   satellites whose ephemerides are in ephem and whose
%                   pseudoranges are in pseudoR
%
%                     [ svID gpsTime elevation azimuth ;
%                       svID gpsTime elevation azimuth ;
%                                    ...
%                       svID gpsTime elevation azimuth ]
%
%         'SVsused' Nsats-by-1 vector of the SV numbers of the
%                   actual satellites that have been used to generate
%                   the navigation solution.  These are the satellites
%                   that obey the elevation mask constraint so that
%                   if mm = find(SVsused(jj,1) == el_az(:,1)),
%                   then el_az(mm,3) >= elevmask.  On the other
%                   hand, if el_az(nn,1) is not an element of
%                   SVsused(:,1), it will be because 
%                   el_az(nn,3) < elevmask at some point of the
%                   solution procedure and probably because
%                   el_az(nn,3) < elevmask on output from this
%                   function.
%
%         'sigmaPR' The RMS pseudorange residual error (meters).
%                   This is the square root of the quantity
%                   derived by dividing the sum of the squares 
%                   of the final pseudorange residuals by 
%                   (Nsats - 4).  This makes sense only if
%                   Nsats > 4.  If Nsats = 4, then this
%                   is output as NaN (not a number) in order to
%                   indicate that there are too few pseudoranges
%                   to compute statistically sensible residuals.
%          
%
function [posOBS,DOP,el_az,SVsused,sigmaPR, ccdot, range, rowindexvec_SVs_used, satX_all, satY_all, satZ_all] = ...
                   solveposod(ephem,pseudoR,guess,gpsTime,...
                              ionParam,iflagion,elevmask,p,...
                              TdegK,hrel,iflagna)
% define physical constants
constant;
% get the number of satellites to use in the navigation solution
Nsatsall = size(ephem,1);
% stop with an error if there are too few satellites.
if Nsatsall < 4
   error(['Error in solveposod.m: fewer than 4',...
          ' satellite pseudoranges entered.'])
end
% compute the transmitter clock times associated with the receiver
% clock time in gpsTime and with the pseudoranges in pseudoR.
% also, ensure that the pseudoranges are stored in a column vector.
pseudoR_colvec_all = pseudoR(:);
t_all = gpsTime - pseudoR_colvec_all*(1/c);
% compute the ECEF positions of the Nsatsall satellites at their 
% transmission times along with their true transmission times and 
% their transmitter clock errors.
satLocClock_all = findsatclock(ephem,t_all);

%grab the time rate of change of the transmitter
%clock error
ccdot=satLocClock_all(:,7);
%%create a row index that will be used for
%%finding the correct beat carrier phase later
Nsats = Nsatsall;
rowindexvec_SVs_used = (1:Nsatsall)';




% compute the time difference vector deltatR_all such that
% deltprop = deltatR_all - recCO is the vector of signal propagation
% delays, in seconds, if recCO is the receiver clock error.
deltatR_all = gpsTime - satLocClock_all(:,2);
% compute pseudoranges that have been corrected for the transmitter
% clock errors.
pseudoCorr_all = satLocClock_all(:,6)*c;
pseudoR_corrected_all = pseudoR_colvec_all + pseudoCorr_all;
% set up the GPS satellite ECEF X and Y positions in the
% transmission time coordinate frames.
satX_Trans_all = satLocClock_all(:,3);
satY_Trans_all = satLocClock_all(:,4);
% set up the GPS satellite ECEF Z positions.
satZ_all = satLocClock_all(:,5);
% initialize arrays that only keep the data for the actual
% satellites that are being used to compute the navigation
% solution, i.e., those whose elevations do not fall
% below the mask limit elevmask.  Initially this consists
% of all the satellite data that have been input.  After
% the solution converges, the elevations will be tested
% against the mask, and satellites will be deleted if
% their elevations fall below the mask.  In that case,
% additional solution iterations are used to refine
% the solution using the new reduced set of satellites.
Nsats = Nsatsall;
SVsused = ephem(:,1);
t = t_all;
deltatR = deltatR_all;
pseudoR_corrected = pseudoR_corrected_all;
satX_Trans = satX_Trans_all;
satY_Trans = satY_Trans_all;
satZ = satZ_all;
% assign the atmospheric parameters to to variable names
% that are used in the call of the neutral atmosphere
% delay model. assign nominal values if values have not  
% been entered. 
if nargin >= 8
   if size(p,1) == 1
      p_used = p;
   else
      p_used = 1013.25;
   end
else
   p_used = 1013.25;
end
if nargin >= 9
   if size(TdegK,1) == 1
      TdegK_used = TdegK;
   else
      TdegK_used = 288.15;
   end
else
   TdegK_used = 288.15;
end
if nargin >= 10
   if size(hrel,1) == 1
      hrel_used = hrel;
   else
      hrel_used = 0.5;
   end
else
   hrel_used = 0.5;
end
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
  % compute the receiver latitude, longitude, and altitude.
  latlongalt = latlong(obsPos');
  % compute the satellite elevations and azimuths.
  el_az = elevazim([SVsused,t,satX,satY,satZ],obsPos');
  elevvec = el_az(:,3);
  % compute the ionospheric delays if they are to be used
  if iflagion == 1
     deltl1ionovec = deltl1ionocalc(el_az,(gpsTime - recCO),...
                                    latlongalt,ionParam);
  else
     deltl1ionovec = zeros(Nsats,1);
  end
  % compute the range-equivalent neutral atmosphere delays if they
  % are to be used.
  if iflagna == 1
     c_deltnavec = deltnacalc(elevvec,p_used,TdegK_used,...
                              hrel_used,latlongalt);
  else
     c_deltnavec = zeros(Nsats,1);
  end
  % form the vector 'l' and the matrix 'A'. Note that the 'A' matrix
  % does not include the small perturbations that would be
  % present if the derivatives of deltl1ionovec with respect to
  % the solved-for user values of X, Y, Z, and recCO and the 
  % derivatives of c_deltnavec with respect to X, Y, and Z.
  % this is reasonable because these perturbations are too small
  % to much affect A or the calculated Newton-Raphson increment x.
  l = pseudoR_corrected - range + arecCO*recCO - c*deltl1ionovec ...
      - c_deltnavec;
  A = [ax,ay,az,(ones(Nsats,1) + arecCO*(1/c))];  
  % solve for "x", which contains deltaX, deltaY, deltaZ, and 
  % the new c*recCO.  Note that this line is equivalent to
  % x = inv(A'*A)*(A'*l) = (A'*A)\(A'*l), that is, left matrix
  % division by a matrix with more rows than columns uses
  % the pseudo-inverse that solves the corresponding least-
  % squares problem.  The actual calculations, however, are
  % somewhat different from inv(A'*A)*(A'*l), different in a
  % way that makes them more stable numerically in the face
  % of the possible poor conditioning of A'*A if GDOP is
  % relatively large.
  x = A\l;
  % Update 'obsPos' by adding the first three elements of 'x' 
  % to the current guess
  obsPos = obsPos + x(1:3,1);
  % Update the receiver clock error estimate.  Also, compute
  % the change in this quantity multiplied by the speed of light.
  delta_c_recCO = x(4,1) - c_recCO;
  c_recCO = x(4,1);
  recCO = c_recCO/c;
  % check to see if the changes are small enough.  if they are,
  % then check whether the elevations of the used satellites
  % are all no less than the mask angle.  if they are,
  % then stop the iteration by setting the 'stop' flag to 1.
  % if the changes in the solution are not small enough, then
  % iterate again.  if the changes are small enough but
  % some of the elevations are too small, then drop the
  % offending satellites and iterate again.
  normerror = norm([x(1:3,1);delta_c_recCO]);
  if normerror < 1.e-06 %1.e-06 m change limit to stop
    ielevdropvec = find(elevvec < elevmask);
    nelevdropvec = size(ielevdropvec,1);
    if nelevdropvec > 0
      Nsats = Nsats - nelevdropvec;
      % stop with an error if there are too few satellites with
      % sufficient elevation.
      if Nsats < 4
         error(['Error in solveposod.m: fewer than 4',...
                ' satellite pseudoranges left after',...
                ' dropping those below the ',num2str(elevmask),...
                ' deg elevation mask angle.'])
      end
      SVsused(ielevdropvec,:) = [];
      t(ielevdropvec,:) = [];
      deltatR(ielevdropvec,:) = [];
      pseudoR_corrected(ielevdropvec,:) = [];
      satX_Trans(ielevdropvec,:) = [];
      satY_Trans(ielevdropvec,:) = [];
      satZ(ielevdropvec,:) = [];
    else
      stop = 1;
    end
  end
  % terminate abnormally if more than 15 iterations are
  % required. send a warning to the display in this case.
  if iters > 15
    stop = 1;
    disp('Warning in solveposod.m: Terminating after 15')
    disp(' iterations without passing convergence test.')
  end
  iters;
end
% create the output matrix 'posOBS'
posOBS = [gpsTime,obsPos',recCO];
% re-compute the satellite azimuths and elevations.
deltheta_all = OmegaE*(deltatR_all - recCO);
cos_deltheta_all = cos(deltheta_all);
sin_deltheta_all = sin(deltheta_all);
satX_all = cos_deltheta_all.*satX_Trans_all + ...
           sin_deltheta_all.*satY_Trans_all;
satY_all = -sin_deltheta_all.*satX_Trans_all + ...
           cos_deltheta_all.*satY_Trans_all;
el_az = elevazim([ephem(:,1),t_all,satX_all,satY_all,satZ_all],obsPos');
% calculate the cofactor matrix 'Q' which is the inverse of the normal
% equation matrix the 'Q' matrix has the following components
% [ qXX qXY qXZ qXt; qYX qYY qYZ qYt; qZX qZY qZZ qZt; qtX qtY qtZ qtt]
Q = inv(A'*A);
% assign diagonal elements 'qXX', 'qYY', 'qZZ', 'qtt'
qXX = Q(1,1);
qYY = Q(2,2);
qZZ = Q(3,3);
qtt = Q(4,4);
% compute 'GDOP', 'PDOP', and 'TDOP' 
GDOP = sqrt(qXX + qYY + qZZ + qtt);
PDOP = sqrt(qXX + qYY + qZZ);
TDOP = sqrt(qtt);
% to compute 'HDOP' and 'VDOP' need rotation matrix from ECEF to local frame
% convert ECEF OBS into latitude-longitude coordinates
psi = latlongalt(1).*degrad;       % latitude
lambda = latlongalt(2).*degrad;    % longitude
% rotation matrix  'R'
R = [ (cos(lambda)*cos(psi)),  (sin(lambda)*cos(psi)), sin(psi); ...
                -sin(lambda),             cos(lambda),        0; ...
     (-cos(lambda)*sin(psi)), (-sin(lambda)*sin(psi)), cos(psi)];
% calculate the local cofactor matrix
Qlocal = [R,zeros(3,1);zeros(1,3),1] * Q * [R',zeros(3,1);zeros(1,3),1];
% assign diagonal elements
qVV = Qlocal(1,1);
qEE = Qlocal(2,2);
qNN = Qlocal(3,3);
% calculate 'HDOP' and 'VDOP' 
HDOP = sqrt(qEE + qNN);
VDOP = sqrt(qVV);	
% return 'DOP'
DOP = [ GDOP; PDOP; TDOP; HDOP; VDOP];
% compute the RMS pseudorange error
if Nsats > 4
  pRresidualsvec = l - A*x;
  sigmaPR = sqrt(sum(pRresidualsvec.^2)/(Nsats - 4));
else
  sigmaPR = NaN;
end
return;