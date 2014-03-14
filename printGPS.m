% Print out nav soln, sat ranges, sat positions, DGPS nav soln, or DOP of a
% nav soln.
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% input: 'arg1' and 'arg2' are the first and second inputs to each print
%        function. The printNav function will run if only one input
%        argument is specified. See individual functions for more details
%        about the input arguments.
%
%        'print' is a string that identifies which print function to run.
%            print = 'ranges' calls printranges
%                  = 'sat'    calls printsat
%                  = 'DGPS'   calls printDGPS
%                  = 'DOP'    calls printDOP

function printGPS(arg1,arg2,print)

if nargin == 1
  printNav(arg1);
elseif strcmp(print,'ranges')
  printranges(arg1,arg2);
elseif strcmp(print,'sat')
  printsat(arg1,arg2);
elseif strcmp(print,'DGPS')
  printDGPS(arg1,arg2);
elseif strcmp(print,'DOP')
  printDOP(arg1,arg2);
elseif strcmp(print,'DOPc')
  outDOPc(arg1,arg2);
else
  fprintf('\nInvaild print function. \n\n')
end

% Print out the navigational solution in ECEF coordinates and
% latitude-longitude-altitude coordinates
%
% input: 'posOBS' matrix which contains a GPS time (seconds), ECEF
%		coordiates of the navigational solution (meters) and the
%		receiver clock offset at that GPS time (seconds)
%						[ GPS time ECEFx ECEFy ECEFz recCO ]
%
function printNav(posOBS)
% get the navigational solution
gpsTime = posOBS(1);
obsECEFx = posOBS(2) / 1000;		% kilometers
obsECEFy = posOBS(3) / 1000; 		% kilometers
obsECEFz = posOBS(4) / 1000;		% kilometers
dT = posOBS(5);
posOBS = latlong(posOBS(2:4));
obsLat = posOBS(1);
obsLong = posOBS(2);
obsAlt = posOBS(3);
% display the navigation solution
fprintf('\n');
fprintf('\nNavigation Solution :');
fprintf('\n*********************\n');
fprintf('\nGPS Time :		  %11.3f  (seconds)',gpsTime);
fprintf('\nLatitude :			%+.5f  (degrees)',obsLat);
fprintf('\nLongitude:			%+.5f  (degrees)',obsLong);
fprintf('\nAltitude :			%+.2f  (meters)',obsAlt);
fprintf('\nClock Offset :			%+.4e  (seconds)',dT);
fprintf('\nX ECEF Coordinate :		%+.3f  (kilometers)',obsECEFx);
fprintf('\nY ECEF Coordinate :		%+.3f  (kilometers)',obsECEFy);
fprintf('\nZ ECEF Coordinate :		%+.3f  (kilometers)',obsECEFz);
% return
fprintf('\n\n');
return


% Print out the pseudo-ranges and calculated ranges to a set of satellites
% over a number of samples
%
function printranges(pseudo,range)
% determine total number of samples
samples = size(pseudo,1);
% determine number of satellites being used
satellites = (size(pseudo,2) - 1) / 2;
% determine which pseudo-range / range samples will be printed
fprintf('\nEnter the pseudo-range/calculated range samples, ');
fprintf('\nto be printed, the samples range from 1 to %d\n', samples);
s = input('	in the form "[ start end ]" or "sample #"  : ');
% determine number of samples to be outputted
num = size(s,2);
if (num == 0)
  s = [ 1 samples ];
else
  if (num == 1)
    samples = 1;
    s(2) = s(1);
  else
    samples = s(2) - s(1) + 1;
  end
end
% print out header
fprintf('\nPseudo-Range / Range sample number');
if (samples == 1)
  fprintf(' %d. ',s(1));
else
  fprintf('s %d to %d. ',s(1),s(2));
end
% print column headers
fprintf('\n\n  GPS Time   ');
for i = 1:satellites
  fprintf('SV ID   Pseudo-Range/Range  ');
end
fprintf('\n');
% print out pseudo-range measurements and range calculations for all
% samples; alternate lines of pseudo-range and range for each sample
for i = s(1):s(2)
  fprintf('%11.3f ',pseudo (i,1));
  for j = 1:satellites
    fprintf('%6d  ',pseudo(i,2 * j ));
    fprintf('%+19.3f ',pseudo(i,2 * j + 1)./1000);
  end
  fprintf('\n            ');
  for j = 1:satellites
    fprintf('%+27.3f ',range(i,2 * j + 1)./1000);
  end
  fprintf('\n');
end
% return
fprintf('\n\n');
return

% Print the calculated satellite positions and elevation-azimuth
%
% input: 'satLoc' matrix which rows contain an SV id number, a GPS
%		time (seconds), and the ECEF coordinates (meters) of the
%		satellite location
%							[ svID GPStime ECEFx ECEFy ECEFz ;
%		  					  svID GPStime ECEFx ECEFy ECEFz ;
%											...
%		  					  svID GPStime ECEFx EFECy ECEFz ]
%        'el-az' martix which rows contain an SV id number, a GPS
%		time (seconds), and the elevation-azimuth look angles
%		(degrees) to the satellite location
%							[ svID GPStime elevation azimuth ;
%		  					  svID GPStime elevation azimuth ;
%											...
%		  					  svID GPStime elevation azimuth ]
%
function printsat(satLoc,el_az)
% find out the time the satellite locations are calculated for; if
% all times are equal, assign 'gpsTime' that value, otherwise assign
% 'gpsTime' equal to zero
gpsTime = satLoc(:,2);
t = sum(gpsTime - gpsTime(1));
if (t == 0)
  gpsTime = gpsTime(1);
else
  gpsTime = 0;
end
% print header to the screen
% check if all satellite positions were found at the same time
fprintf('\n\nSatellite locations ');
if (gpsTime ~= 0)
  fprintf('for GPS time of week  : %11.3f sec',gpsTime);
end
% print column headers
fprintf('\n\n%4s  %11s  %11s  %11s  %8s  %8s  %11s\n', ...
  'SV', 'X (km)', 'Y (km)', 'Z (km)', 'el (deg)', ...
  'az (deg)', 'r (km)');
% output satellite's elevation-azimuth and ECEF coordinates
satX = satLoc(:,3) ./ 1000;		% convert to kilometers
satY = satLoc(:,4) ./ 1000;		% convert to kilometers
satZ = satLoc(:,5) ./ 1000;		% convert to kilometers
range = sqrt((satX.^2) + (satY.^2) + (satZ.^2)); %this is from center of earth...
for i = 1:size(satLoc,1)
  % print SV number
  fprintf('%4d  ',satLoc(i,1));
  % print ECEF coordinates in kilometers
  fprintf('%+11.3f  %+11.3f  %+11.3f  ', ...
    satX(i),satY(i),satZ(i));
  % print elevation and azimuth
  fprintf('%+8.1f  %+8.1f  ',el_az(i,3),el_az(i,4));
  % print satellite range
  fprintf('%11.3f\n', range(i));
end
% return
fprintf('\n\n')
return

% Prints out the navigational solution in latitude-longitude-altitude
% coordinates for both uncorrected and DGPS  corrected pseudo-ranges
%
function printDGPS(posOBS,posOBSwDGPS)
% get the navigational solution of the uncorrected pseudo-range solution
gpsTime = posOBS(1);
posOBS = latlong(posOBS(2:4));
obsLat = posOBS(1);
obsLong = posOBS(2);
obsAlt = posOBS(3);
% display the navigation solution for uncorrected pseudo-range solutions
fprintf('\n');
fprintf('\nNavigation Solution with Uncorrected Pseudo-Ranges:');
fprintf('\n***********************************************\n');
fprintf('\nGPS Time :         %11.3f  (seconds)',gpsTime);
fprintf('\nLatitude  :        %+.5f  (degrees)',obsLat);
fprintf('\nLongitude  :       %+.5f  (degrees)',obsLong);
fprintf('\nAltitude  :        %+.2f  (meters)',obsAlt);
% get the navigational solution of the DGPS corrected solution
gpsTime = posOBSwDGPS(1);
posOBSwDGPS = latlong(posOBSwDGPS(2:4));
obsLat = posOBSwDGPS(1);
obsLong = posOBSwDGPS(2);
obsAlt = posOBSwDGPS(3);
% display the navigation solution for the DGPS corrected solutions
fprintf('\n');
fprintf('\nNavigation Solution with DGPS Corrected Pseudo-Ranges:');
fprintf('\n***************************************************\n');
fprintf('\nGPS Time :         %11.3f  (seconds)',gpsTime);
fprintf('\nLatitude  :        %+.5f  (degrees)',obsLat);
fprintf('\nLongitude  :       %+.5f  (degrees)',obsLong);
fprintf('\nAltitude  :        %+.2f  (meters)',obsAlt);
% return
fprintf('\n\n');
return

% Print out the dilution of precision of a navigational solution
function printDOP(posOBS,DOP)
% get the navigational solution
pos=latlong(posOBS(:,2:4));
gpsTime = posOBS(1);
obsLat = pos(1);
obsLong = pos(2);
obsAlt = pos(3);
% get the dilution of precision calculations
GDOP = DOP(1);
PDOP = DOP(2);
TDOP = DOP(3);
HDOP = DOP(4);
VDOP = DOP(5);
% display the navigation solution
fprintf('\n');
fprintf('\nNavigation Solution :');
fprintf('\n*********************\n');
fprintf('\nGPS Time  :		    %11.3f (seconds)',gpsTime);
fprintf('\nLatitude  :          %+.5f  (degrees)',obsLat);
fprintf('\nLongitude :          %+.5f  (degrees)',obsLong);
fprintf('\nAltitude  :			 %+.2f    (meters)',obsAlt);
fprintf('\n\nDilution of Precision :');
fprintf('\n*********************\n');
fprintf('\nGeomteric Dilution of Precision :	 %.3f',GDOP);
fprintf('\nPositional Dilution of Precision :	 %.3f',PDOP);
fprintf('\nTime Dilution of Precision :		 %.3f',TDOP);
fprintf('\nHorizontal Dilution of Precision :	 %.3f',HDOP);
fprintf('\nVertical Dilution of Precision :	 %.3f',VDOP);
fprintf('\n\n');
return

% outDOPc.m
%
% Prints the results of navDOP in tabular form including
% Combination Number, SVIDs, Error, PDOP and summary statistics.
% Also plots a comparision of PDOP and Error
%
% input: 'obsArray' containing the navigation solution, error, PDOP
%				etc. See navDOP for a full description
%			'combination' a vector containing all possible 4
%				satellite combinations for a given sample number
%
function outDOPc(obsArray,combination)
% print 4 satellite results
num=(1:size(combination,1))';
fprintf('\n');
fprintf('\n4 Satellite Solution Results:');
fprintf('\n*****************************\n');
fprintf('\nCombo#  SVID1  SVID2  SVID3  SVID4  Error(meters)  PDOP');
for i=1:size(num,1)
  fprintf('\n  %2.0f', num(i,1))
  fprintf('      %2.0f', combination(i,1))
  fprintf('     %2.0f', combination(i,2))
  fprintf('     %2.0f', combination(i,3))
  fprintf('     %2.0f', combination(i,4))
  fprintf('       %6.1f', obsArray(i,6))
  fprintf('    %5.1f', obsArray(i,5))
end;
fprintf('\n\n');
% calculate and print 4-satellite solution accuracy statistics
meanError=mean(obsArray(:,6));
stdError=std(obsArray(:,6));
minError=min(obsArray(:,6));
maxError=max(obsArray(:,6));
fprintf('\nMaximum Error :   %6.1f meters :  4 SVIDs', maxError);
fprintf('\nMinimum Error :   %6.1f meters :  4 SVIDs', minError);
fprintf('\nMean Error    :   %6.1f meters :  4 SVIDs', meanError);
fprintf('\n\n');
% plot position error and PDOP vs. combination number
figure
subplot(2,1,1);
stem(1:size(combination,1),obsArray(:,6));
title('4 Satellite Navsoln Error/PDOP vs. SV Combination Number')
ylabel('Error (m)'); xlim([0,size(combination,1)+1])
subplot(2,1,2);
stem(1:size(combination,1),obsArray(:,5));
ylabel('PDOP'); xlim([0,size(combination,1)+1])
xlabel('Combination #');
return