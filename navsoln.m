% NAVSOLN.m		(actual file name: NAVSOLN.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% < FOUR SATELLITE SOLUTION >
%
% This program determines a navigational solution directly from raw
% pseudo-ranges and satellite ephemerides. The user is required to 
% supply an initial receiver location (within 100 kilometers), four
% satellite SV ids to use in the navigational solution from, and a
% pseudo-range measurement sample which specifies a GPS time
% and pseudo-range measurements to be used to calculate the
% navigational solution.
%
% The files required for NAVSOLN.m are:
%		* constant.m
%		* ecef.m
%		* findsatclock.m
%		* formatdata.m
%		* latlong.m
%		* printGPS.m
%		* solvepos.m
%		
% Ephemeris data and pseudo-range data are needed for calculating
% the corrected pseudo-ranges, the satellite positions, and the
% final navigational solution. The program parseRinex.m will format
% data into two corresponding files; ephem.asc and obs.asc. 
%
% NAVSOLN.m consists of the following steps and function calls:		
%		* opens data files - ephem.asc and obs.asc
%		* inputs the observation station location
%		* inputs satellite four SVs to be used
%		* calls formatdata to format 'ephem' and 'obs' data
%			structures
%		* inputs which pseudo-range measurement is to be used
%			to calculate the position
%		* determine pseudo-range measurements and GPS time of
%			sample chosen
%		* calls solvepos which iteratively calculates the receiver's
%			position from satellite locations and pseudo-ranges 	
%		* calls printGPS to print out results 
%
% < FOUR SATELLITE SOLUTION >
%
% retrieve ephemeris data from input file -- ephem.asc 
	load ephem.asc;
	ephemData = ephem; 
	clear ephem;
% retrieve pseudo-range data from input file -- pseudo.asc
	load obs.asc; 
	pseudoData = obs;
	clear obs; 
% input initial guess of receiver's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate location of the observation station \n');
	guess = input('	in the form "[ latitude longitude altitude ]"  :  ');
    guess = ecef(guess);
% get the list of SVIDs for which obs and ephem data are available.
    SVIDlist = pseudoData(:,3:2:end);
    SVIDlist = unique(SVIDlist(:));
    if SVIDlist(1,1) == 0
       SVIDlist(1,:) = [];
    end
    nsatsdum = size(SVIDlist,1);
    for k = nsatsdum:-1:1
       if ~any(SVIDlist(k,1) == ephemData(:,1))
          SVIDlist(k,:) = [];
       end
    end
    clear nsatsdum;
% determine which satellites will be used in the navigational
% solution
	fprintf('\nEnter the satellites to be used in the navigational solution ');
	fprintf('\nEnter FOUR, no less, no more - they must have pseudoranges:');
	fprintf('\nfrom the list of SV ids --> ');
    fprintf('%d ',SVIDlist);
    fprintf('\n');
	SV = input('	in the form of   "[ sv sv . . . ]"  :  ');
% call formatData which will reformat 'ephemData' and 'pseudoData'
% into the structures 'ephem' and 'pseudo' respectively
	[ ephem pseudo ] = formatdata(ephemData,pseudoData,SV);
    pseudo
% check that four satellites have been specified
	if (size(ephem,1) ~= 4)
		fprintf('\nNeed to specify four satellites with ephemerides !!!\n\n');
		return;
	end
% determine which pseudo-range measurement will be used
	fprintf('\nEnter the pseudo-range measurement sample, which range');
	fprintf('\nfrom 1 to %d, to be used to calculate the position \n',size(pseudo,1));
	s = input('	in the form of  "sample #"    :  ');
% check that valid sample number has been entered
	if ((s < 1) || (s > size(pseudo,1)))
		fprintf('\nSample number is out of range !!!\n\n');
		return;
	end
% determine pseudo-range measurements 'pseudoR' and GPS time 
% 'gpsTime' of the sample chosen
	pseudoR = [ pseudo(s,3) ; pseudo(s,5); ...
			    pseudo(s,7); pseudo(s,9) ];
	gpsTime = pseudo(s,1);
% check that have psuedo-range data for all satellites for the 
% sample
	if (~isempty(find(pseudoR == 0,1)))
		fprintf('\nIncomplete pseudo-range data for the sample !!!\n\n');
		return;
	end
% call solvePos passing pseudo-ranges 'pseudoR', an initial 
% positional guess 'guess', a GPS time 'gpsTime', and ephemeris 
% data 'ephem'
	posOBS = solvepos(ephem,pseudoR,guess,gpsTime);
% output the navigational solution
	printGPS(posOBS);		    