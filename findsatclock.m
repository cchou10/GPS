% findsatclock.m	(actual file name: findsatclock.m)
%
%  Copyright (c) 2008 Paul M. Kinter, Jr., and 
%  Copyright (c) 2012 Mark L. Psiaki.  All rights reserved. 
%
%	< INCLUDES KEPLERIAN CORRECTIONS AND DETERMINES SATELLITE
%     CLOCK ERROR >
%
% this GPS function computes satellite positions given the orbital
% ephemerides at given times expressed as (erroneous) transmitter
% clock times.  It also compute the corresponding true times
% along with the the transmitter clock errors at those times.
%
% input: 'ephem' satellites-by-24 matrix, each of whose rows contains 
%                orbital ephemerides and other navigation message
%                data for a given satellite.

%					< see formatData.m for description >

%	 		't'  satellites-by-1 vector, each of whose rows contains  
%                the GPS times, according to the corresponding
%                satellite's transmitter clock, at which the true 
%                transmission time, the position, and the satellite 
%                clock error will be calculated.  Note: unlike
%                findsat.m, this function requires that the number
%                of rows in t and the number of rows in ephem
%                be identical.  Also unlike findsat.m, this time
%                input is an erroneous transmitter clock time
%                rather than a true time.
%
% output: 'satLocClock'  satellites-by-6 matrix each of whose rows contain
%                        an SV id number, the true GPS time (seconds)
%                        that corresponds to the GPS transmitter clock 
%                        time in the corresponding row of t, the 
%                        ECEF coordinates (meters) of the satellite 
%		                 location in the ECEF coordinate system 
%                        that applies at the true transmission time, 
%                        and the transmitter clock error (seconds)

%						 [ svID GPStime ECEFx ECEFy ECEFz cc;
%		  				   svID GPStime ECEFx ECEFy ECEFz cc;
%											...
%						   svID GPStime ECEFx EFECy ECEFz cc]
%
%                        Note: the relationship between the true
%                        times in satLocClock(:,2), the satellite
%                        transmitter clock errors in satLocClock(:,6),
%                        and the satellite clock times in t is the 
%                        following:
%
%                           satLocClock(:,2) = t - satLocClock(:,6)
%
%	< INCLUDES KEPLERIAN CORRECTIONS AND DETERMINES SATELLITE
%     CLOCK ERROR >
%
function satLocClock = findsatclock(ephem,t)
% define physical constants
constant;
% determine number of satellites; return an error message if zero satellites
satellites = size(ephem,1);
if (satellites == 0)
   error('Number or rows in ephem is zero in findsatclock.m.')
end
% ensure that the number of rows in t equals the number of rows in ephem
if size(t,1) ~= satellites
   error('Number or rows in t and in ephem do not agree in findsatclock.m.')
end
% define orbital parameters
t0e = ephem(:,4);		% ephemeris reference time (seconds)
ecc = ephem(:,5);		% eccentricity (unitless)
sqrta = ephem(:,6);		% square root of semi-major axis (meters^(1/2))
omega0 = ephem(:,7);	% argument of perigee (radians)
M0 = ephem(:,8);		% mean anomaly at reference time (radians)
l0 = ephem(:,9);		% right ascension at reference (radians)
lDot = ephem(:,10);     % rate of right acension (radians/second)
dn = ephem(:,11);		% mean motion difference (radians/second)
i0 = ephem(:,12); 		% inclination angle at reference time (radians)
iDot = ephem(:,13); 	% inclination angle rate (radians/second)
cuc = ephem(:,14);		% latitude cosine harmonic correction (radians)
cus = ephem(:,15);		% latitude sine harmonic correction (radians)
crc = ephem(:,16);		% orbit radius cosine harmonic correction (meters)
crs = ephem(:,17);		% orbit radius sine harmonic correction (meters)
cic = ephem(:,18);		% inclination cosine harmonic correction (radians)
cis = ephem(:,19);		% inclination sine harmonic correction (radians)
% get clock correction parameters from 'ephem'
refTime = ephem(:,24);  % clock correction reference time (seconds)
af0 = ephem(:,20);      % clock correction bias term (seconds)
af1 = ephem(:,21);      % clock correction rate term (seconds/sec)
af2 = ephem(:,22);      % clock correction acceleration term (seconds/sec^2)
tgd = ephem(:,23);      % single-frequency additional clock corr. (seconds)
% compute constant coefficient of relativistic clock correction term
Frel_sqrta_ecc = (- 2*sqrt(muearth)/(c^2))*(sqrta.*ecc);
% develop first guess of true transmission time in a two-step
% process
T = t - af0 + tgd;
deltaTweek = -604800*round((T - refTime)*(1/604800));
T = (t - af0 - af1.*(deltaTweek - refTime) + tgd)./(1 + af1);
% develop first guess of eccentric anomaly as mean anomalay
% at first guess of true transmission time.
n0 = sqrt(muearth) * ((sqrta).^(-3));
n = n0 + dn;
T_corr = T - 604800*round((T - t0e)*(1/604800));
dt = T_corr - t0e;
M = M0 + n.*dt;
E = M;
% compute the true transmission time and the eccentric anomaly 
% by simultaneously solving the equations:
%		f1(T,E) = T - t + cc(T,E)
%       f2(T,E) = M(T) - E + ecc*sin(E)
% where cc(T,E) is the transmitter clock error
for k = 1:10
    % compute f1 and f2 values
    deltaTweek = -604800*round((T - refTime)*(1/604800));
    timeOffset = T + deltaTweek - refTime;
    sin_E = sin(E);
    cos_E = cos(E);
    cc = af0 + af1.*timeOffset + af2.*(timeOffset.^2) - tgd ...
           + Frel_sqrta_ecc.*sin_E;
    f1 = T - t + cc;    
    T_corr = T - 604800*round((T - t0e)*(1/604800));
    dt = T_corr - t0e;
    M = M0 + n.*dt;
    f2 = M - E + ecc.*sin_E;
    % compute partial derivatives of f1 and f2 with respect to T and E.
    dccdT = af1 + 2*(af2.*timeOffset);
    df1dT = 1 + dccdT;
    dccdE = Frel_sqrta_ecc.*cos_E;
    df1dE = dccdE;
    df2dT = n;
    df2dE = - 1 + ecc.*cos_E;
    % solve for Newton-Raphson increment using a method
    % that is valid for each separate satellite's linear
    % system of 2 equations in 2 unknowns.
    determinantvec = df1dT.*df2dE - df1dE.*df2dT;
    oodeterminantvec = determinantvec.^(-1);
    dT = - (df2dE.*f1 - df1dE.*f2).*oodeterminantvec;
    dE = - (- df2dT.*f1 + df1dT.*f2).*oodeterminantvec;
    % perform Newton-Raphson updates.
    T = T + dT;
    E = E + dE;
    % end early if all eccentric anomaly changes produce less
    % than about 0.001 m of satellite position change and if
    % all clock corrections produce less than 0.001 m of satellite
    % position error at a velocity of 4000 m/sec.
    if (max(abs(dE)) <= 3.5e-11) & (max(abs(dT)) <= 2.5e-07)
        break 
    end
end
% re-calculate quantities used in T and E iterative calculation
% one final time to be sure they are at their final converged values.
deltaTweek = -604800*round((T - refTime)*(1/604800));
timeOffset = T + deltaTweek - refTime;
sin_E = sin(E);
cos_E = cos(E);
cc = af0 + af1.*timeOffset + af2.*(timeOffset.^2) - tgd ...
           + Frel_sqrta_ecc.*sin_E;
T_corr = T - 604800*round((T - t0e)*(1/604800));
dt = T_corr - t0e;

% calculate true anomoly from eccentric anomoly
sinnu = (sqrt(1 - ecc.^2) .* sin_E) ./ (1 - ecc .* cos_E);
cosnu = (cos_E - ecc) ./ (1 - ecc .* cos_E);
nu = atan2(sinnu,cosnu);
% calculate the argument of latitude
Phi = nu + omega0;
% calculate the corrected argument of perigee.
twoPhi = 2*Phi;
cos2Phi = cos(twoPhi);
sin2Phi = sin(twoPhi);
omegaCorr = cuc .* cos2Phi + cus .* sin2Phi;
omega = omega0 + omegaCorr;
% calculate longitude of ascending node with correction
lcorr = lDot.*dt;
l = l0 - OmegaE*T_corr + lcorr;
% calculate orbital radius with correction
rCorr = crc .* cos2Phi + crs .* sin2Phi;
r = (sqrta.^2) .* (1 - ecc .* cos_E) + rCorr;
% calculate inclination with correction
iCorr = cic .* cos2Phi + cis .* sin2Phi + iDot .* dt;
i = i0 + iCorr;
% find position in orbital plane
u = omega + nu;
xp = r .* cos(u);
yp = r .* sin(u);
% find satellite position in ECEF coordinates
cos_i = cos(i);
sin_i = sin(i);
cos_l = cos(l);
sin_l = sin(l);
yp_cos_i = yp .* cos_i;
ECEFx = (xp .* cos_l) - (yp_cos_i .* sin_l);
ECEFy = (xp .* sin_l) + (yp_cos_i .* cos_l);
ECEFz = (yp .* sin_i);


%%calculate the transmitter clock error time derivative
%%first calculating Mdot then Edot
Mdot = n;
Edot = Mdot./(1 - ecc.*cos_E);
% compute the transmitter clock error time derivatives
ccdot = af1 + 2*af2.*timeOffset + Frel_sqrta_ecc.*cos_E.*Edot;





% return true transmission times, satellite locations, and transmitter
% clock errors.
satLocClock = [ ephem(:,1), T, ECEFx, ECEFy, ECEFz, cc, ccdot];
return;