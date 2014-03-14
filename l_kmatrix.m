% l_kmatrix.m	(actual file name: l_kmatrix.m)
%
% this CDGPS utility function computes the l_k matrix at a particular
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
%           'pos_known'             position known in ECEF coordinates of
%                                   the known receiver
%
% output:   'l_k'                   the l_k matrix used for various
%                                   caculations. size: [N-1] by 1
%
%           'range_known'           the distance from the receiver to the 
%                                   satellite based only on receiver and
%                                   satellite position
%
function [l_k range_known] = l_kmatrix(pseudoRA, pseudoRB,...
                    ephem, BCphase_A, BCphase_B,gpsTimeA,gpsTimeB, obsPos,...
                    ionParam,iflagion,elevmask,p,...
                              TdegK,hrel,iflagna, pos_known)
% call solveposod_n for data set A 
    [posOBS_A,DOP_A,el_az_A,SVsused_A,sigmaPR_A, ccdotA, range_A, rowindexvec_SVs_used, satX, satY, satZ] = ...
        solveposod_n(ephem,pseudoRA,obsPos,gpsTimeA,...
                              ionParam,1,elevmask,...
                              p,TdegK,hrel,1);
% call solveposod_n for data set B
    [posOBS_B,DOP_B,el_az_B,SVsused_B,sigmaPR_B, ccdotB, range_B, rowindexvec_SVs_used, satX_known, satY_known, satZ_known] = ...
        solveposod_n(ephem,pseudoRB,pos_known,gpsTimeB,...
                              ionParam,1,elevmask,...
                              p,TdegK,hrel,1);                         
% calculate the true range from the known reciever based
% on the known position of the reciever station and the satellites
% positions at the time of transmission
    deltaXsatrcvr_known = satX_known - pos_known(1);
    deltaYsatrcvr_known = satY_known - pos_known(2);
    deltaZsatrcvr_known = satZ_known - pos_known(3);
    range_known = sqrt(deltaXsatrcvr_known.^2 + deltaYsatrcvr_known.^2 + deltaZsatrcvr_known.^2);
% load the constant fule
    constant;
% get the true time of reception for data set A and B
    t_receptionA = posOBS_A(:,1)-posOBS_A(:,5);
    t_receptionB = posOBS_B(:,1)-posOBS_B(:,5);                            
% define the difference in reception time
    delta_t_reception = t_receptionA-t_receptionB;
% define the satellites beat carrier phases
% that will actually be used for calculating
% the relative navigational solution
    BCphase_A = BCphase_A(rowindexvec_SVs_used,1);
    BCphase_B = BCphase_B(rowindexvec_SVs_used,1);
% define single difference beat carrier phase
    delta_BC_A = BCphase_A(1)-BCphase_A(2:end);
    delta_BC_B = BCphase_B(1)-BCphase_B(2:end);
% define double difference beat carrier phase
    delta_delta_BC = delta_BC_A-delta_BC_B;
% compute the ranges to the satellites.
    deltaXsatrcvr = satX - obsPos(1);
    deltaYsatrcvr = satY - obsPos(2);
    deltaZsatrcvr = satZ - obsPos(3);
    range = sqrt(deltaXsatrcvr.^2 + deltaYsatrcvr.^2 + deltaZsatrcvr.^2);
% define the coefficients
    freq_co1 = 1+ccdotA(1);
    freq_coN = 1+ccdotA(2:end);
% define the range from reciever A and B to the first satellite
    rangeA_1 = range(1);
    rangeB_1 = range_B(1);
% define the ranges from reciever's A and B to 2-N satellites
    rangeA_2N = range(2:end);
    rangeB_2N = range_B(2:end);
% define the rate of change of transmitter clock error
% from satellite 1 to reciever A and
% from satellites 2-N to reciever A
    ccdotA1 = ccdotA(1);
    ccdotA2N = ccdotA(2:end);
    l_k = lambdaL1.*delta_delta_BC - freq_co1.*(rangeA_1 - rangeB_1)...
        + freq_coN.*(rangeA_2N - rangeB_2N)...
      +c*(ccdotA1 - ccdotA2N)*delta_t_reception;  
end