% poserror.m	(actual file name: poserror.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
% calculates the error in the navigation solution by finding the
% vector distance between the calculated and actual positions
%
% input: 'calculated' calculated ECEF navigation solution
%			'acutal' actual ECEF coordinated of the antenna
%					[ calculated, actual ]
%
function distance = poserror(calculated, actual)
% calculate the difference
xDiff = calculated(1)-actual(1);
yDiff = calculated(2)-actual(2);
zDiff = calculated(3)-actual(3);
% square the difference
xSqr = xDiff^2;
ySqr = yDiff^2;
zSqr = zDiff^2;
squares = [xSqr,ySqr,zSqr];
% calculate the error
distance = sqrt(sum(squares));
return