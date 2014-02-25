% deltnacalc.m	(actual file name: deltnacalc.m)
%
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved.  
%
%	< COMPUTES NEUTRAL ATMOSPHERE DELAY CORRECTIONS >
%
% this function calculates the pseudo-range corrections
% that compensate for the delay caused by the Earth's
% neutral atmosphere.  It uses the Saastamoinen zenith
% dry and wet delays and the Ifadis dry and wet elevation
% mapping functions.
%
% inputs:   'Elvec'  The Nsats-by-1 vector of elevations of the
%                    line-of-sight vectors to the GPS satellites
%                    whose neutral atmosphere delays are to
%                    be calculated (degrees)
%
%                      [ El; El; El; ... ]
%
%		    'p'      Atmospheric pressure at receiver (millibar or,
%                    identically, hPa)
%
%           'TdegK'  Temperature at receiver (degrees Kelvin)
%
%           'hrel'   Relative humidity at receiver (fraction in 
%                    range 0 to 1)
%
%           'ecoord' 1-by-3 vector which contains the position specified
%                    by WGS-84 latitude (degrees), longitude (degrees),  
%                    and altitude (meters)
%
%                      [ latitude, longitude, altitude ]
%
% outputs:  'c_deltnavec' The Nsats-by-1 vector of neutral atmosphere
%                         delays multiplied by the speed of light 
%                         (meters) along the line-of-sight paths
%                         whose elevations are given in Elvec.
%
%                           [ c_deltna; c_deltna; c_deltna; ... ]
%
function c_deltnavec = deltnacalc(Elvec,p,TdegK,hrel,ecoord)
% Compute the water vapor partial pressure in hPa units.
% Use pp. 36-38, 41, 42 of the Mendes dissertation.
   TdegK0 = 273.15;         % degrees K
   esat0 = 6.11;            % mb or hPa
   Lv = 2.50e+06;           % Joules/kg
%  Ridealgas = 8314.510;    % Joules/(kmole-degK)
%  Mwater = 18.01528;       % kg/kmole
%  Rw = Ridealgas/Mwater;
   Rw = 461.525;            % Joules/(kg-degK)
   ooTdegK = 1/TdegK;
   esat = esat0*exp((Lv/Rw)*((1/TdegK0) - ooTdegK));
   TdegC = TdegK - 273.15;
   fw = 1.00072 + p*(3.2e-06 + 5.9e-10*(TdegC^2));
   esatpr = esat*fw;
   e = esatpr*hrel/(1 - (1 - hrel)*(esatpr/p));
% Compute the range-equlvalent Saastamoinen zenith dry/hydrostatic 
% delay as a function of pressure, latitutde, and altitude, in meters.
   degrad = pi/180.0;
   twolat = 2*ecoord(1,1)*degrad;
   costwolat = cos(twolat);
   deltaroh_zd = 0.002277*p/(1 - 0.0026*costwolat - 0.00000028*ecoord(1,3));
% Compute the range-equlvalent Saastamoinen zenith non-hydrostatic 
% delay (the wet delay) as a function of temperature and water 
% vapor partial pressure, in meters.
   deltaroh_zw = 0.002277*(1255*ooTdegK + 0.05)*e;
% Compute some quantities that will be useful for computing
% the hydrostatic and non-hydrostatic mapping functions.
   pm1000 = p - 1000;
   TdegKm288 = TdegK - 288.15;
   sqrte = sqrt(e);
% Compute the sines of the elevation angles, but keep the elevations
% from going to zero or below by clipping the elevations
% at 0.1 deg if they are less than 0.1 deg.  This clipping
% prevents division by zero in the mapping function calculations.
% it may be very helpful when calling this function from within
% a navigation solution algorithm that may start with a very
% bad first guess of position.
   Elvec_clipped = Elvec;
   ilowvec = find(Elvec_clipped < 0.1);
   nlow = size(ilowvec,1);
   if nlow > 0
      Elvec_clipped(ilowvec,1) = 0.1*ones(nlow,1);
   end
   sinElvec = sin(Elvec_clipped*degrad);
% Compute the coefficients for the Ifadis hydrostatic elevation 
% mapping function.
   ad = 0.001237 + (0.1316e-06)*pm1000 + (0.1378e-05)*TdegKm288 + ...
                    (0.8057e-05)*sqrte;
   bd = 0.003333 + (0.1946e-06)*pm1000 + (0.1040e-06)*TdegKm288 + ...
                    (0.1747e-04)*sqrte;
   cd = 0.078;
% Compute the Ifadis hydrostatic elevation mapping function at
% the various elevations in Elvec.
   numfac1d = 1 + cd;
   oonumfac1d = 1/numfac1d;
   numfac2d = 1 + bd*oonumfac1d;
   oonumfac2d = 1/numfac2d;
   numfac3d = 1 + ad*oonumfac2d;
   denfac1dvec = sinElvec + cd;
   oodenfac1dvec = denfac1dvec.^(-1);
   denfac2dvec = sinElvec + bd*oodenfac1dvec;
   oodenfac2dvec = denfac2dvec.^(-1);
   denfac3dvec = sinElvec + ad*oodenfac2dvec;
   oodenfac3dvec = denfac3dvec.^(-1);
   mapdvec = numfac3d*oodenfac3dvec;
% Compute the coefficients for the Ifadis non-hydrostatic 
% (i.e., wet) elevation mapping function.
   aw = 0.0005236 + (0.2471e-06)*pm1000 - (0.1724e-06)*TdegKm288 + ...
                     (0.1328e-04)*sqrte;
   bw = 0.001705 + (0.7384e-06)*pm1000 + (0.3767e-06)*TdegKm288 + ...
                    (0.2147e-04)*sqrte;
   cw = 0.05917;
% Compute the Ifadis non-hydrostatic (i.e., wet) elevation mapping 
% function at the various elevations in Elvec
   numfac1w = 1 + cw;
   oonumfac1w = 1/numfac1w;
   numfac2w = 1 + bw*oonumfac1w;
   oonumfac2w = 1/numfac2w;
   numfac3w = 1 + aw*oonumfac2w;
   denfac1wvec = sinElvec + cw;
   oodenfac1wvec = denfac1wvec.^(-1);
   denfac2wvec = sinElvec + bw*oodenfac1wvec;
   oodenfac2wvec = denfac2wvec.^(-1);
   denfac3wvec = sinElvec + aw*oodenfac2wvec;
   oodenfac3wvec = denfac3wvec.^(-1);
   mapwvec = numfac3w*oodenfac3wvec;
%
%  Finish computing the total range-equivalent of the
%  neutral atmosphere delay by multiplying the mapping functions
%  by the zenith delays and by adding the hydrostatic and
%  non-hydrostatic (wet) components together.
% 
   c_deltnavec = deltaroh_zd*mapdvec + deltaroh_zw*mapwvec;