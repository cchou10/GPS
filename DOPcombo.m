% DOPcombo.m		(actual file name: DOPcombo.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% This program examines the relationship between navigation solution
% accuracy and DOP.
%
% The files required for DOPcombo.m are:
%	all files for 'NAVSOLN.m' plus
%		* satPos.m
%		* combo.m
%		* choose.m
%		* printGPS.m
%		* calcdop.m
%		* poserror.m
%

% clear the Matlab workspace
clear
% retrieve ephemeris data from input file -- ephem.asc
load ephem.asc;
ephemData = ephem;
clear ephem;

% retrieve pseudo-range data from input file -- pseudo.asc
load obs.asc;
pseudoData = obs;
clear obs;
% get the receiver antenna location from the user
guess = input('Enter antenna location in the form [lat (deg) long (deg) alt (m)]:\n');
guess = ecef(guess);
% determine which pseudo-range measurement will be used
fprintf('\nEnter the pseudo-range measurement sample, which range');
fprintf('\nfrom 1 to %d, to be used to calculate the position ',size(pseudoData,1));
s = input('\n	in the form of  "sample #"    :  ');
% check that valid sample number has been entered
if ((s < 1) || (s > size(pseudoData,1)))
  fprintf('\nSample number is out of range\n\n');
  return;
end
% call formatData which will reformat 'ephemData' and 'pseudoData'
% into the structures 'ephem' and 'pseudo' respectively
[ ephem pseudo ] = formatData(ephemData,pseudoData,ephemData(:,1)');
% call pseudoCalc passing 'ephem', 'pseudo', 'ionParam', and 'guess'
% which contains the raw pseudo-ranges, satellite clock corrections,
% ionospheric parameters, and observation location guess
pseudo_range = pseudoCalc(ephem,pseudo(s,:));
% determine pseudo-range measurements 'pseudoR' for all available
% satellites and GPS time 'gpsTime' of the sample chosen
SVID_list = pseudoData(s,3:2:end);
SVID_list = unique(SVID_list(:));
if SVID_list(1,1) == 0
   SVID_list(1,:) = [];
end
nsatsdum = size(SVID_list,1);
for k = nsatsdum:-1:1
  if ~any(SVID_list(k,1) == ephem(:,1))
    SVID_list(k,:) = [];
  end
end
% get the raw pseudorange data that go with this list of SVs in
% SVID_list.
nsatsdum = size(SVID_list,1);
pseudo_list = zeros(nsatsdum,1);
pseudoData_SVs = pseudoData(s,3:2:end)';
pseudoData_PRs = pseudoData(s,4:2:end)';
for k = 1:nsatsdum
   lldumk = find(SVID_list(k,1) == pseudoData_SVs);
   pseudo_list(k,1) = pseudoData_PRs(lldumk,1);
end
clear nsatsdum pseudoData_SVs pseudoData_PRs lldumk
gpsTime = pseudo_range(1);
% if less than 4 satellites, stop calculation
if length(SVID_list)<4
  error('Less than 4 satellites at this sample number');
end
% print out satellite position
satPos(ephem,gpsTime,guess,SVID_list)
% determine which satellites will be used in the navigational
% solution
fprintf('\nEnter the satellites to be used in the navigational solution(s) ');
fprintf('\nfrom the list of SV ids --> ');
fprintf('%d ',SVID_list);
fprintf('\n');
SV = input('	in the form of   "[ sv sv . . . ]"  :  ');
% default is all satellites
if isempty(SV)
  SV = SVID_list;
end
% check that four satellites have been specified
if (length(SV) < 4)
  fprintf('\nNeed to specify at least four satellites\n\n');
  return
end
% get all possible combinations of four satellites
combination=combo(SV,4);
% calculate navigation solution, positional error, pdop and gdop
% for all combinations of 4 satellites
obsArray=[];
posOBSmat = [];
for i=1:size(combination,1)
  ephemCombo=[];
  for j=1:size(ephem,1)
    if (~isempty(find(combination(i,1:4)==ephem(j,1), 1)))
      ephemCombo=[ephemCombo;ephem(j,:)];
    end
  end
  pseudoCombo=zeros(4,1);
  for j=1:4
    lldumj = find(ephemCombo(j,1) == SVID_list);
    pseudoCombo(j,1) = pseudo_list(lldumj,1);
  end
  clear lldumj
  % find calculated user position for the current combination of SVs
  posOBS = solvepos(ephemCombo,pseudoCombo,guess,gpsTime);
  posOBSmat = [posOBSmat;posOBS];
  % calculate DOP for the current combination of SVs
  iflagA = 0;
  DOP = calcdop(ephemCombo,gpsTime-posOBS(1,5),posOBS(1,2:4)',iflagA);
  % calculate error for each navitgation solution
  positionError = poserror(posOBS(2:4),guess);
  % save the results in obsArray
  obsArray=[obsArray;posOBS(2),posOBS(3),posOBS(4),...
    DOP(1),DOP(2),positionError];
end
% print results
printGPS(obsArray,combination,'DOPc')