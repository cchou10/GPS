
function Qtilda=calculateQtilda(pseudoRA, pseudoRB,...
                    ephem, BCphase_A, BCphase_B,gpsTime, guess,...
                    ionParam,iflagion,elevmask,p,...
                              TdegK,hrel,iflagna)
guessB=guess;
guessA=guess;
gpsTimeA=gpsTime;
gpsTimeB=gpsTime;
ephemA=ephem;
ephemB=ephem;
 %call solve posod for data set A        
                
[posOBS_A,DOP_A,el_az_A,SVsused_A,sigmaPR_A, ccdotA, range_A, rowindexvec_SVs_used] = ...
    solveposod(ephemA,pseudoRA,guessA,gpsTimeA,...
                              ionParam,1,elevmask,...
                              677,1,1,1);
                  
                
                
  %call solve posod for data set B
[posOBS_B,DOP_B,el_az_B,SVsused_B,sigmaPR_B, ccdotB, range_B, rowindexvec_SVs_used] = ...
                   solveposod(ephemB,pseudoRB,guessB,gpsTimeB,...
                              ionParam,1,elevmask,...
                              667,1,1,1);      
                
%%calculate the true time of transmission
%%for reciever A
cc_dotA=ccdotA;


%%get the true time of reception for data set A and B
t_receptionA=posOBS_A(:,1)-posOBS_A(:,3);
t_receptionB=posOBS_B(:,1)-posOBS_B(:,3);                
                
                
%%define the difference in reception time
delta_t_reception=t_receptionA-t_receptionB;

%%define the satellites beat carrier phases
%%that will actually be used for calculating
%%the relative navigational solution
BCphase_A_used=BCphase_A(rowindexvec_SVs_used,1);
BCphase_B_used=BCphase_B(rowindexvec_SVs_used,1);

% define physical constants
constant;
% compute the transmitter clock times associated with the receiver
% clock time in gpsTime and with the pseudoranges in pseudoR.
% also, ensure that the pseudoranges are stored in a column vector.
pseudoR_colvec = pseudoRA(:);
t = gpsTime - pseudoR_colvec*(1/c);
% compute the ECEF positions of the 4 satellites at their transmission
% times along with their true transmission times and their transmitter
% clock errors.
satLocClock = findsatclock(ephem,t);
% compute the time difference vector deltatR such that
% deltprop = deltatR - recCO is the vector of signal propagation
% delays, in seconds, if recCO is the receiver clock error.
deltatR = gpsTime - satLocClock(:,2);
% compute pseudoranges that have been corrected for the transmitter
% clock errors.
pseudoCorr = satLocClock(:,6)*c;
pseudoR_corrected = pseudoR_colvec + pseudoCorr;
% set up the GPS satellite ECEF X and Y positions in the
% transmission time coordinate frames.
satX_Trans = satLocClock(:,3);
satY_Trans = satLocClock(:,4);
% set up the GPS satellite ECEF Z positions.
satZ = satLocClock(:,5);
% initialize the guesses of the receiver position and of the
% receiver clock error.  also initialize a variable that holds
% the receiver clock error multiplied by the speed of light.
obsPos = guess(:);
recCO = min(deltatR) - 0.067;
c_recCO = c*recCO;
% initialize iteration counter.
iters=0;
% solve for position iteratively until solution is within an
% acceptable error
%stop = 0;
%while (stop == 0)
 % iters=iters+1;
 
 
  % rotate satellite position vectors into ECEF reference frame
  % of the current guess of the true reception time.
  deltheta = OmegaE*(deltatR - recCO);
  cos_deltheta = cos(deltheta);
  sin_deltheta = sin(deltheta);
  satX = cos_deltheta.*satX_Trans + sin_deltheta.*satY_Trans;
  satY = -sin_deltheta.*satX_Trans + cos_deltheta.*satY_Trans;
  % compute the ranges to the satellites.
  deltaXsatrcvr = satX - obsPos(1,1);
  deltaYsatrcvr = satY - obsPos(2,1);
  deltaZsatrcvr = satZ - obsPos(3,1);
  range = sqrt(deltaXsatrcvr.^2 + deltaYsatrcvr.^2 + deltaZsatrcvr.^2);
   % compute the partial derivatives of the ranges with respect to the
  % elements of obsPos.
  oorange = range.^(-1);
  ax = - deltaXsatrcvr.*oorange;
  ay = - deltaYsatrcvr.*oorange;
  az = - deltaZsatrcvr.*oorange;
  
  %define the partial derivatives with respect to 
  %one satellite then with respect to the rest of
  %the satellites
  ax_1=ax(1);
  ay_1=ay(1);
  az_1=az(1);
  
  ax_N=ax(2:end);
  ay_N=ay(2:end);
  az_N=az(2:end);
  
  
  %define the coefficients
  freq_co1=1+cc_dotA(1);
  freq_coN=1+cc_dotA(2:end);
  
  %define single difference beat carrier phase
  delta_BC_A=BCphase_A(1)-BCphase_A(2:end);
  delta_BC_B=BCphase_B(1)-BCphase_B(2:end);
  %define double difference beat carrier phase
  delta_delta_BC=delta_BC_A-delta_BC_B;
  
  
  %define the range from reciever A and B to the first satellite
  rangeA_1=range(1);
  rangeB_1=range_B(1);
  %define the ranges from reciever's A and B to 2-N satellites
  rangeA_2N=range(2:end);
  rangeB_2N=range_B(2:end);
  %define the rate of change of transmitter clock error
  %from satellite 1 to reciever A
  ccdotA1=cc_dotA(1);
  
  %define the rate of change of the transmitter clock error
  %from satellites 2-N to reciever A
  ccdotA2N=cc_dotA(2:end);
  
  %%define the A matrix
  A=[freq_co1.*ax_1-freq_coN.*ax_N...
     freq_co1.*ay_1-freq_coN.*ay_N...
     freq_co1.*az_1-freq_coN.*az_N];
  
  
 

  %%define the l matrix

  l=lambdaL1.*delta_delta_BC...
      -freq_co1.*(rangeA_1-rangeB_1)...
      +freq_coN.*(rangeA_2N-rangeB_2N)...
      +c*(ccdotA1-ccdotA2N)*delta_t_reception;
  
  %%define the Q matrix
  
  Q_dd=2.*ones(length(A))+2.*eye(length(A));
  
  %%define the qtilda matrix which will be the output and used
  %to calculate the biases
  Qtilda=inv(inv(Q_dd)-inv(Q_dd)*A*inv(A'*inv(Q_dd)*A)*A'*inv(Q_dd));
  

  
  
  
end