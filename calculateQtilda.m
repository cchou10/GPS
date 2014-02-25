% calculateQtilda.m	(actual file name: calculateQtilda.m)
%
% this CDGPS utility function computes inverse Qtilda at a particular
% sample
%
% input:    'pseudoRA'              pseudorange from obsA.asc for receiver A formatted so
%                                   its N by 1 with each element the 
%                                   pseudorange to a satellite 
%
%           'pseudoRB'              pseudorange from obsB.asc for receiver B formatted so
%                                   its N by 1 with each element the 
%                                   pseudorange to a satellite.
%
%           'ephem'                 formatted ephemeride data for each
%                                   satellite
%
%           'BCphase_A'             the beat carrier phases for each
%                                   satellite from receiver A
%
%           'BCphase_B'             beat carrier phases for each satellite
%                                   from receiver B
%
%           'gpsTimeA'              time observed by receiver A. number of
%                                   sesconds since Sunday midnight
%
%           'gpsTimeB'              time observed by receiver B. number of
%                                   sesconds since Sunday midnight
%
%           'obsPos'                observed position of the receivers in
%                                   ECEF form
%
%           'ionoParam'             the ionosphere measurements loaded from
%                                   ion.asc
%
%           'iflagion'              1 if including ionosphere calculations
%                                   and 0 if not
%
%           'elevmask'              threshold of which satellites to use in
%                                   degrees
%
%           'p'                     environmental pressure at that time
%                                   outside in mmHg
%
%           'TdegK'                 degrees celsius of the temperature 
%
%           'hrel'                  the humidity between 0 and 1 of the
%                                   environment
%
%           'iflagna'               1 if including neutral atmosphere calculations
%                                   and 0 if not
%
% output:   'Qtilda_inv'            the calculated value of Qtilda_inv that
%                                   is used to calculate the real-value bias 
%                                   and Jsql
%
%           'A'                     A_k matrix that is used to for calculation 
%                                   of Qtilda as well as Jsql, biases, and
%                                   delta X, Y, and Z
%
function [Qtilda_inv A]=calculateQtilda(pseudoRA, pseudoRB,...
                    ephem, BCphase_A, BCphase_B,gpsTime, obsPos,...
                    ionParam,iflagion,elevmask,p,...
                              TdegK,hrel,iflagna)
% calls solveposod         
    [posOBS_A,DOP_A,el_az_A,SVsused_A,sigmaPR_A, ccdotA, range_A, rowindexvec_SVs_used,satX, satY, satZ] = ...
                    solveposod_n(ephem,pseudoRA,obsPos,gpsTime,...
                              ionParam,1,elevmask,...
                              677,1,1,1);     
% compute the ranges to the satellites.
    deltaXsatrcvr = satX - obsPos(1);
    deltaYsatrcvr = satY - obsPos(2);
    deltaZsatrcvr = satZ - obsPos(3);
    range = sqrt(deltaXsatrcvr.^2 + deltaYsatrcvr.^2 + deltaZsatrcvr.^2);
% compute the partial derivatives of the ranges with respect to the
% elements of obsPos.
    oorange = range.^(-1);
    ax = - deltaXsatrcvr.*oorange;
    ay = - deltaYsatrcvr.*oorange;
    az = - deltaZsatrcvr.*oorange;
% define the partial derivatives with respect to 
% one satellite then with respect to the rest of
% the satellites
    ax_1 = ax(1);
    ay_1 = ay(1);
    az_1 = az(1);  
    ax_N = ax(2:end);
    ay_N = ay(2:end);
    az_N = az(2:end); 
% define the coefficients
    freq_co1 = 1+ccdotA(1);
    freq_coN = 1+ccdotA(2:end);  
% define the A matrix
    A = [freq_co1.*ax_1-freq_coN.*ax_N...
        freq_co1.*ay_1-freq_coN.*ay_N...
        freq_co1.*az_1-freq_coN.*az_N];
% define the Q matrix
    Q_dd = 2.*ones(length(A))+2.*eye(length(A));
% define the qtilda matrix which will be the output and used
% to calculate the biases
    invQdd = inv(Q_dd);
    Qtilda_inv = invQdd-(invQdd*A*inv(A'*invQdd*A)*A'*invQdd);  
end