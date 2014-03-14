% calcdop.m		(actual file name: calcdop.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved.
%
% this function calculates the dilution of precision of a 
% navigational solution based on the satellite geometry, given the
% navigation data, the true reception time, and the
% receiver true position.
%
% input: 'ephem'  Nsats-by-24 matrix, each of whose rows contains 
%                 orbital ephemerides and other navigation message
%                 data for a given satellite.
%
%					< see formatdata.m for description >
%
%		 'gpsTime' variable which contains the true GPS time of
%                  reception (seconds) for which the dilution of
%                  precision is being calculated.
%
%		 'obsPos'  3-by-1 vector containing the known position of
%                  the receiver in ECEF Cartesian coordinates (meters)
%                  at which the dilution of precision is to
%                  be calculated.
%
%        'iflagA'  The scalar input flag that determines whether
%                  the dilution of precision calculations use
%                  the standard A matrix or whether the last column
%                  of the A matrix includes the small effects of
%                  receiver clock time uncertainty on the
%                  propagation delay calculation in the 
%                  position/clock solution and, therefore, on
%                  the dilution of precison.  If iflagA = 0, then
%                  the standard A matrix is used, but if
%                  iflagA = 1, then the modified A matrix
%                  is used, the one that includes 
%
% output: 'DOP'    5-by-1 column vector which contains the various 
%                  measures of dilution of precision
%
%                     DOP(1,1) geometrical dilution of precision
%                     DOP(2,1) positional dilution of precision
%                     DOP(3,1) time dilution of precision
%                     DOP(4,1) horizontal dilution of precision
%                     DOP(5,1) vertical dilution of precision
% 
function DOP = calcdop(ephem,gpsTime,obsPos,iflagA)
  constant;
% determine number of satellites
  Nsats = size(ephem,1);
% initialize the propagation delays.
  deltprop = zeros(Nsats,1);
% iterate to calculate the satellite positions along with the
% delay propagation.  this is an ad hoc iteration that converges
% very rapidly due to the fact that the speed of light is larger
% than the inertial velocities of the GPS satellites by about
% 5 orders of magnitude.
  for k = 1:7
% calculate the satellite locations at times of transmission
% in ECEF frames at times of transmission
    satXYZ = findsat(ephem,(gpsTime - deltprop));
	satX_Trans = satXYZ(:,3);
	satY_Trans = satXYZ(:,4);
	satZ = satXYZ(:,5);
% rotation to ECEF frame at time of reception.
    deltheta = OmegaE*deltprop;
    cos_deltheta = cos(deltheta);
    sin_deltheta = sin(deltheta);
    satX = cos_deltheta.*satX_Trans + sin_deltheta.*satY_Trans;
    satY = -sin_deltheta.*satX_Trans + cos_deltheta.*satY_Trans;
% compute the ranges to the satellites.
    deltaXsatrcvr = satX - obsPos(1);
    deltaYsatrcvr = satY - obsPos(2);
    deltaZsatrcvr = satZ - obsPos(3);
    range = sqrt(deltaXsatrcvr.^2 + deltaYsatrcvr.^2 + deltaZsatrcvr.^2);
% calculate the new propagation times and quit early if they
% differ from the old ones by less than 1.e-07 sec.
    deltprop_old = deltprop;
    deltprop = range*(1/c);
    if max(abs(deltprop - deltprop_old)) < 1.e-07
       break
    end
  end
% compute the partial derivatives of the ranges with respect to the
% elements of obsPos.
  oorange = range.^(-1);
  ax = - deltaXsatrcvr.*oorange;
  ay = - deltaYsatrcvr.*oorange;
  az = - deltaZsatrcvr.*oorange;
% compute the A matrix.  branch depending on how the propagation 
% delay is to be treated in computing the last column of the A matrix.
  if iflagA == 0
    A = [ax,ay,az,ones(Nsats,1)];  
  else
% compute the partial derivatives of the ranges with respect recCO.
    arecCO = OmegaE*((-deltaXsatrcvr.*satY + deltaYsatrcvr.*satX).*oorange);
%  form the vector 'l' and the matrix 'A'
    A = [ax,ay,az,(ones(Nsats,1) + arecCO*(1/c))];  
  end
% calculate the cofactor matrix 'Q' which is the inverse of the normal
% equation matrix the 'Q' matrix has the following components
% [ qXX qXY qXZ qXt; qYX qYY qYZ qYt; qZX qZY qZZ qZt; qtX qtY qtZ qtt]
  Q = inv(A' * A);
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
  obsPos_row = obsPos(:);
  obsPos_row = obsPos_row';
  latlongalt = latlong(obsPos_row);
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
return