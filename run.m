% run.m		(actual file name: run.m)
%
% This program calculates things for double differencing 
%
% The files required are:
%
%		* constant.m
%       * deltl1ionocalc.m
%       * deltnacalc.m
%		* ecef.m
%       * elevazim.m
%		* findsatclock.m
%		* formatdata.m
%		* latlong.m
%		* printGPS.m
%		* solveposod.m
%		
% Ephemeris data, pseudo-range data, and ionosphere model data
% are needed for calculating the corrected pseudo-ranges, the 
% satellite positions, and the final navigational solution. 
% The program parseRinex.m must have already formated data into  
% three corresponding files; ephem.asc, ion.asc, and obs.asc. 
%
% run.m consists of the following steps and function calls:
%
%		* opens data files - ephem.asc, ion.asc, and obs.asc
%       * inputs neutral atmosphere pressure (p), temparature (TdegC),
%         and relative humidity (hrel) data and computes temperature
%         in degrees Kelvin (TdegK).  If empty arrays are input, 
%         then nominal values will be used: p = 1013.25 millibars,
%         TdegK =  288.15 deg K, and hrel = 0.50.
%       * inputs flags that determine whether ionosphere delay
%         corrections and neutral atmosphere delay corrections
%         are applied, during both solution for clock errors
%         and calculation of pseudorange errors.
%		* inputs the observation station location guess
%       * inputs elevation mask angle.
%		* inputs four or more satellite SVs to be used
%		* calls formatdata to format 'ephem' and 'obs' data
%			structures
%		* inputs which pseudo-range measurement is to be used
%			to calculate the position
%		* determine pseudo-range measurements and GPS time of
%			sample chosen
%		* calls solveposod which iteratively calculates the receiver's
%			position from satellite locations and pseudo-ranges 
%       * prints out satellite usage, elevation, and azimuth information
%         along with DOP and pseudorange RMS error information
%		* calls printGPS to print out results 
%
%

% clear the Matlab workspace
    clear
% retrieve ephemeris data from input file -- ephemA.asc for first receiver
% data
	load ephem.asc;
	ephemData = ephem; 
	clear ephem;
% retrieve ionosphere model data from input file -- ion.asc for first
% receiver data
    load ion.asc
    ionParam = ion;
    clear ion
% retrieve pseudo-range data from input file -- pseudo.asc for first
% receiver data
	load obsA.asc; 
	pseudoDataA = obsA;
	clear obsA;  
% retrieve pseudo-range data from input file -- pseudo.asc for second
% receiver data
	load obsB.asc; 
	pseudoDataB = obsB;
	clear obsB;
% Load the beat carrier for receiver A
    load('L1GPS3.mat');
    phiA_data = obs;
    clear obs;
% Load the beat carrier phase for receiver B
    load('L1GPS6.mat');
    phiB_data = obs;
    clear obs;
% input elevation mask angle. 
	fprintf('\nEnter the elevation mask angle below which ');
	fprintf('\n the correspoonding pseudorange will not \n');
	elevmask = input(' be used (deg) elevmask :  ');
% input initial guess of receiver A's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate location of the observation station A \n');
	guess = input(' in the form "[ latitude longitude altitude ]"  :  ');
    guessA = ecef(guess);
% get the list of SVIDs for which obs and ephem data are available.
    SVIDlistA = pseudoDataA(:,3:2:end);
    SVIDlistA = unique(SVIDlistA(:));
    if SVIDlistA(1,1) == 0
       SVIDlistA(1,:) = [];
    end
    nsatsdum = size(SVIDlistA,1);
    for k = nsatsdum:-1:1
       if ~any(SVIDlistA(k,1) == ephemData(:,1))
          SVIDlistA(k,:) = [];
       end
    end
    clear nsatsdum;
% determine which satellites will be used in the navigational
% solution
	fprintf('\nEnter the satellites to be used in the navigational solution ');
	fprintf('\nEnter FOUR or more - they must have pseudoranges:');
	fprintf('\nfrom the list of SV ids --> ');
    fprintf('%d ',SVIDlistA);
    fprintf('\n');
	SVA = input('	in the form of   "[ sv sv . . . ]"  :  ');
% input initial guess of receiver's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate location of the observation station B \n');
	guess = input(' in the form "[ latitude longitude altitude ]"  :  ');
    guessB = ecef(guess);
% get the list of SVIDs for which obs and ephem data are available.
    SVIDlistB = pseudoDataB(:,3:2:end);
    SVIDlistB = unique(SVIDlistB(:));
    if SVIDlistB(1,1) == 0
       SVIDlistB(1,:) = [];
    end
    nsatsdum = size(SVIDlistB,1);
    for k = nsatsdum:-1:1
       if ~any(SVIDlistB(k,1) == ephemData(:,1))
          SVIDlistB(k,:) = [];
       end
    end
    clear nsatsdum;
% determine which satellites will be used in the navigational
% solution
	fprintf('\nEnter the satellites to be used in the navigational solution ');
	fprintf('\nEnter FOUR or more - they must have pseudoranges:');
	fprintf('\nfrom the list of SV ids --> ');
    fprintf('%d ',SVIDlistB);
    fprintf('\n');
	SVB = input('	in the form of   "[ sv sv . . . ]"  :  ');
% call formatData which will reformat 'ephemData' and 'pseudoData'
% into the structures 'ephem' and 'pseudo' respectively
	[ ephemA pseudoA ] = formatdata(ephemData,pseudoDataA,SVA);
    [ ephemB pseudoB ] = formatdata(ephemData,pseudoDataB,SVB);
    
%call formatData which will reformat 
    [ ephemA phiA ] = formatdata(ephemData,phiA_data,SVA);
    [ ephemB phiB ] = formatdata(ephemData,phiB_data,SVB);
% check that four satellites have been specified
	if (size(ephemA,1) < 4 || size(ephemB,1) < 4)
		fprintf(['\nNeed to specify four or more satellites',...
                 ' with ephemerides !!!\n\n']);
		return;
    end
% determine which pseudo-range measurement will be used
	fprintf('\nEnter the pseudo-range measurement sample, which range');
	fprintf('\nfrom 1 to %d, to be used to calculate the position \n',...
            size(pseudoA,1));
	sA = input('	in the form of  "sample #"    :  ');
    fprintf('\nEnter the pseudo-range measurement sample, which range');
	fprintf('\nfrom 1 to %d, to be used to calculate the position \n',...
            size(pseudoB,1));
	sB = input('	in the form of  "sample #"    :  ');
% check that valid sample number has been entered
	if ((sA < 1) || (sB <1) || (sA > size(pseudoA,1)) || sB > size(pseudoB,1))
		fprintf('\nSample number is out of range !!!\n\n');
		return;
	end
% determine pseudo-range measurements 'pseudoR' and GPS time 
% 'gpsTime' of the sample chosen
	pseudoRA = pseudoA(sA,3:2:end)';
    pseudoRB = pseudoB(sB,3:2:end)';
%deterime the beat carrier phase
%measurements and GPS time of the choosen sample
    phiA=phiA(sA,3:2:end)';
    phiB=phiB(sB,3:2:end)';
    
% check that have psuedo-range data for all satellites for the 
% sample
    if any(pseudoRA == 0)
       idumbadvec = find(pseudoRA == 0);
       SVbadvec = pseudoA(sA,2:2:end)';
       SVbadvec = SVbadvec(idumbadvec,1);
       disp(' ')
       disp('Error in navsolnod.m: The following SVs')
       disp(' with receiver A appear to lack valid pseudorange data:')
       fprintf('%d ',SVbadvec);
       disp(' ')
       return;
    end
	gpsTimeA = pseudoA(sA,1);
    
    if any(pseudoRB == 0)
       idumbadvec = find(pseudoRB == 0);
       SVbadvec = pseudoB(sB,2:2:end)';
       SVbadvec = SVbadvec(idumbadvec,1);
       disp(' ')
       disp('Error in navsolnod.m: The following SVs')
       disp(' with receiver B appear to lack valid pseudorange data:')
       fprintf('%d ',SVbadvec);
       disp(' ')
       return;
    end
	gpsTimeB = pseudoB(sB,1);
% call solveposed passing pseudo-ranges 'pseudoR', an initial 
% positional guess 'guess', a GPS time 'gpsTime', and ephemeris 
% data 'ephem', and other relevant parameters. This is for receiver A
    [posOBSA,DOPA,el_azA,SVsusedA,sigmaPRA] = ...
                   solveposod(ephemA,pseudoRA,guessA,gpsTimeA,...
                              ionParam,1,elevmask,...
                              677,1,1,1);
% output the navigational solution
	printGPS(posOBSA);
    
% call solveposed passing pseudo-ranges 'pseudoR', an initial 
% positional guess 'guess', a GPS time 'gpsTime', and ephemeris 
% data 'ephem', and other relevant parameters. This is for receiver B
    [posOBSB,DOPB,el_azB,SVsusedB,sigmaPRB] = ...
                   solveposod(ephemB,pseudoRB,guessB,gpsTimeB,...
                              ionParam,1,elevmask,...
                              667,1,1,1);
% output the navigational solution
	printGPS(posOBSB);

    
    

obsPos_noway=solveDDGPS(pseudoRA, pseudoRB,...
                    ephemA, phiA, phiB,gpsTimeA, guess,...
                    ionParam,1,elevmask,113,...
                              677,1,1);
                          
                  ecoord = latlong(obsPos_noway')
