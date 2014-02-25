% satPos.m  (actual file name: satpos.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% plots the location of all tracked satellites at given time
% and prints a table of SVID elevaion and azimuth data
%
% input:	'ephem' from formatda output, 
%			'gpsTime' the GPS time for the calculation
%        'obsLoc' the precise antenna location (ECEF) 
%        'SVID' the list of satellites to be plotted
%						[ ephem, gpsTime, obsLoc, SVID ]
%
% requires: plotsat.m and all associated subprograms 
%
function [] = satPos(ephem,gpsTime,obsLoc,SVID)
% find satellite locations using findsat.m
	satLoc = findsat(ephem,gpsTime);
% build elevation/azimuth table using elevazim.m and sort by SVID
	unsort_matrix = elevazim(satLoc,obsLoc);
	matrix = [];
	for i = 1:32
		for j = 1:(size(unsort_matrix,1))
			if unsort_matrix(j,1) == i
				matrix = [matrix;unsort_matrix(j,:)];
			end
		end
	end
% restrict to tracked satellites
	data=[];
	for i=1:size(matrix,1)
		if ~isempty(find(SVID==matrix(i,1), 1))
			data=[data;matrix(i,:)];
		end
	end	
% print sorted elevation/azimuth values in a table
	fprintf('\n');
	fprintf('\nSatellite Elevation/Azimuth Data:');
	fprintf('\n*********************************\n');
	fprintf('\nSVID  Elevation(deg) Azimuth(deg)');
	for i=1:size(data,1)
		fprintf('\n %2.0f', data(i,1))
		fprintf('       %5.1f', data(i,3))
		fprintf('         %5.1f', data(i,4))
		end;
	fprintf('\n\n');
% plot satellite postions
	plotsat(data)
 	return