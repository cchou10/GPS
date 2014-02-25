% deltl1ionocalc.m	(actual file name: deltl1ionocalc.m)
%
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved.  
%
%	< COMPUTES IONOSPHERE DELAY CORRECTIONS >
%
% this function calculates the modeled ionospheric 
% delays that will be used to compensate the pseudoranges
% for the delays caused by the Earth's ionosphere.  It
% uses the Klobuchar model whose parameters are broadcast
% in the GPS navigation data message.
%
% inputs:   'el_az'     The Nsats-by-4 matrix whose rows contain  
%                       an SV id number, a GPS time (seconds), 
%                       and the elevation and azimuth look angles
%		                (degrees) to the satellite whose L1
%                       ionospheric delay will be returned
%                       in the corresponding row of deltl1ionovec
%
%						   [ svID GSPtime elevation azimuth ;
%						     svID GSPtime elevation azimuth ;
%                                         ...
%		  				     svID GSPtime elevation azimuth ]
%
%           TR          The reception time of week (seconds) at
%                       which the delays are to be computed.  If
%                       possible, this should be the true
%                       GPS time of week.
%
%           'ecoord'    1-by-3 vector which contains the position 
%                       specified by WGS-84 latitude (degrees), 
%                       longitude (degrees), and altitude (meters)
%
%                         [ latitude, longitude, altitude ]
%
%           'ionParam'  4-by-2 matrix of the eight coefficients 
%                       used to correct for ionospheric delay
%                       (seconds in first row, seconds/semi-circle
%                       in second row, seconds/(semi-circle)^2 in
%                       third row, and seconds/(semi-circle)^3 in
%                       4th row.)
%
%                          [ alpha0 beta0
%                            alpha1 beta1
%                            alpha2 beta2
%                            alpha3 beta3 ]
%
% outputs:  'deltl1ionovec' The Nsats-by-1 vector of ionospheric
%                           delays (seconds) at the GPS L1
%                           frequency along the line-of-sight paths
%                           whose elevations and azimuths are
%                           are given in el_az.
%
%                              [ deltl1iono; deltl1iono; deltl1iono; ... ]
%
function deltl1ionovec = deltl1ionocalc(el_az,TR,ecoord,ionParam)
% Compute the great circle arcs from the receiver points to
% the ionospheric pierce points.  before calculating,
% clip the elevation to be no lower than 0 deg.  this
% is a precaution against bad effects from the possibility
% of very negative elevations being used when this function
% is called from within a navigation solution that starts
% with a very poor initial guess of its position.
   degrad = pi/180.0;
   elradians = el_az(:,3)*degrad;
   idum = find(elradians < 0);
   ndum = size(idum,1);
   if ndum > 0
      elradians(idum,1) = zeros(ndum,1);
   end
   psi = 0.1352*((elradians + 0.346).^(-1)) - 0.0691;
% Compute the geodetic latitudes of the ionospheric pierce
% points.
   azradians = el_az(:,4)*degrad;
   phii = ecoord(1,1)*degrad + psi.*cos(azradians);
   ilowvec = find(phii < -1.307);
   nlow = size(ilowvec,1);
   if nlow > 0
      phii(ilowvec,1) = -1.307*ones(nlow,1);
   end
   ihivec = find(phii > 1.307);
   nhi = size(ihivec,1);
   if nhi > 0
      phii(ihivec,1) = 1.307*ones(nhi,1);
   end
% Compute the geodetic longitudes of the ionospheric pierce points.
   lambdai = ecoord(1,2)*degrad + ((psi.*sin(azradians))./cos(phii));
% Compute the geomagnetic latitudes of the ionospheric pierce points.
   phim = phii + 0.201*cos(lambdai - 5.08);
% Translate the geomagnetic latitudes into semi-circle units and
% compute their 2nd and 3rd powers.
   phim_semi = phim*(1/pi);
   phim_semisq = phim_semi.^2;
   phim_semicu = phim_semisq.*phim_semi;
% Compute the amplitude coefficients.
   C1 = ionParam(1,1) + ionParam(2,1)*phim_semi + ...
        ionParam(3,1)*phim_semisq + ionParam(4,1)*phim_semicu;
   ilowvec = find(C1 < 0);
   nlow = size(ilowvec,1);
   if nlow > 0
      C1(ilowvec,1) = zeros(nlow,1);
   end
% Compute the cosinusoidal variation period.
   Per = ionParam(1,2) + ionParam(2,2)*phim_semi + ...
         ionParam(3,2)*phim_semisq + ionParam(4,2)*phim_semicu;
   ilowvec = find(Per < 72000);
   nlow = size(ilowvec,1);
   if nlow > 0
      Per(ilowvec,1) = 72000*ones(nlow,1);
   end
% Compute the local times of day at the pierce points.
   ti = mod((TR + (43200/pi)*lambdai),86400);
% Compute the pierce point longitude angles relative to 2 p.m.
% local time.
   theta = (2*pi)*((ti - 50400)./Per);
% Compute the approximation cos(theta), but zero it out if
% the absolute value of theta is above 1.57.
   costhetaapprox = 1 - 0.5*(theta.^2) + (1/24)*(theta.^4);
   ihivec = find(abs(theta) >= 1.57);
   nhi = size(ihivec,1);
   if nhi > 0
      costhetaapprox(ihivec,1) = zeros(nhi,1);
   end
% Compute the vertical ionospheric delays, i.e., without the slant
% factor.
   v_deltl1ionovec = 5.e-09 + C1.*costhetaapprox;
% Compute the slant factors using the approximate formula
   F = 1 + (16/(pi^3))*((0.53*pi - elradians).^3);
% Use the slant factors and the vertical delays to compute
% the slant delays.
   deltl1ionovec = F.*v_deltl1ionovec;
return