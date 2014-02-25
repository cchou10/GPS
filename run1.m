% run.m		(actual file name: run1.m)
%
% This program calculates the position of a reciever by implementing the
% double-differenced CDGPS algorithm
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
%		* solveposod.m
%       * hyperrectangle.m
%       * calculateQtilda.m
%       * l_kmatrix.m
%       * calculateBiases.m
%       * solveDDGPS.m
%		
% Ephemeris data, pseudo-range data, and ionosphere model data
% are needed for calculating the corrected pseudo-ranges, the 
% satellite positions, and the final navigational solution through CDGPS. 
% The program parseRinex.m must have already formated data into  
% three corresponding files; ephem.asc, ion.asc, and obs.asc. 
%
% run.m consists of the following steps and function calls:
%
%		* opens data files - ephem.asc, ion.asc, obsA.asc, and obsB.asc
%       * opens beat carrier phase data - L1GPS3.mat and L1GPS6.mat
%		* inputs the observation station location guess for both receivers
%		* inputs four or more satellite SVs to be used for both receivers
%		* calls formatdata to format 'ephem' and 'obs' data structures
%		* inputs which pseudo-range measurement is to be used
%			to calculate the position
%		* determine pseudo-range measurements and GPS time of
%			sample chosen
%		* calls calculateQtilda which calculates inverse Qtilda
%           the A_k matrix for each sample time
%       * calls l_kmatrix which calculates the l_k matrix for each sample
%           time
%       * using inverse Qtilda and l_k, calculate a real-value bias estimate
%       * using hyperrectangle.m, calculate the integer-ambiguity lattice
%           points. inputs theoretical standard deviation and aeta are
%           calculated as well
%       * iterate through to determine the best integer-ambiguity lattice
%           points
%       * call solveDDGPS.m to determine the change of position and iterate
%           through until maximum accuracy is achieved
%
%

% clear the Matlab workspace
    clear
    
fprintf('choose which static reciever data set you would like to use\n');
fprintf('1- Rhodes hall North and Rhodes hall south. RHS is static reciever\n');
fprintf('2- Upson Hall and Rhodes hall south. Upson Hall is static reciever\n');
rec_choice=input('input the data set number (1 or 2):');



if(rec_choice==1)
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
	pseudoDataB = obsA;
	clear obsA;  
% retrieve pseudo-range data from input file -- pseudo.asc for second
% receiver data
	load obsB.asc; 
	pseudoDataA = obsB;
	clear obsB;
% load the beat carrier phase for each receiver
    load('L1GPS3.mat');
    phiA_data = obs;
    clear obs;
    load('L1GPS6.mat');
    phiB_data = obs;
    clear obs;
    % input the known reciever location
    pos_known=ecef([42.443327 -76.481433 250.335]);
    %input the surveyed solution to the unknown reciever
    pos_unknown_rcvr=ecef([42.443545 -76.481790 253.94]);
    
elseif(rec_choice==2)
    % retrieve ephemeris data from input file -- ephemA.asc for first receiver
% data
	load pre_ephem.asc;
	ephemData = pre_ephem; 
	clear ephem;
% retrieve ionosphere model data from input file -- ion.asc for first
% receiver data
    load ion_data2.asc
    ionParam = ion_data2;
    clear ion
% retrieve pseudo-range data from input file -- pseudo.asc for first
% receiver data
	load ('gps1Obs.mat'); 
	pseudoDataA = C1;
    phiA_data = L1;
	clear L1;
    clear C1;
% retrieve pseudo-range data from input file -- pseudo.asc for second
% receiver data
	load ('gpsrec01Obs.mat'); 
	pseudoDataB = C1;
    phiB_data = L1;
	clear L1;
    clear C1;
    % input the known reciever location
    pos_known=ecef([42.444007 -76.482229 236.548]); 
    %input the surveyed solution to the unknown reciever
    pos_unknown_rcvr=ecef([42.443327 -76.481433 250.335]);
    
end
    
    
% input elevation mask angle. 
	elevmask = 0;
% input initial guess of receiver A and B's position in latitude-longitude
% coordinates 
	fprintf('\nEnter the approximate location of the unknown observation station \n');
	guess = input(' in the form "[ latitude longitude altitude ]"  :  ');
    guessA = ecef(guess);
    guessB = ecef(guess);
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
% solution for both receivers
if(rec_choice==1)
	SVA = [1 11 12 14 22 25 30 31 32];
    SVB = SVA;
elseif(rec_choice==2)
    SVA = [2 4 5 10 12 13 25];
    SVB = SVA;
end


% call formatData which will reformat 'ephemData' and 'pseudoData'
% into the structures 'ephem' and 'pseudo' respectively
	[ephemA pseudoA] = formatdata(ephemData,pseudoDataA,SVA);
    [ephemB pseudoB] = formatdata(ephemData,pseudoDataB,SVB);
% call formatData which will reformat 'ephemData' and 'phi_data' into the
% structures 'ephemA' and 'phi' respectively
    [ephemA phiA] = formatdata(ephemData,phiA_data,SVA);
    [ephemB phiB] = formatdata(ephemData,phiB_data,SVB);
% check that four satellites have been specified
	if (size(ephemA,1) < 4 || size(ephemB,1) < 4)
		fprintf(['\nNeed to specify four or more satellites',...
                 ' with ephemerides !!!\n\n']);
		return;
    end
    
    
%align the pseudorange and dopplershift data for option 2
    if(rec_choice==2)
    phiA=phiA(98:end,:);
    pseudoA=pseudoA(98:end,:);
    end
% define the sample size to calculate biases and relative position
    n = 2000;
% define a matrix that will be used to store all the current positions
% from the different time iterations
    obsPos_matrix = [];
% define the guess as the current observation position of the known
% receiver
    obs = guessB;
% create a for loop to add in the guesses of all the current positions
    for i = 1:n
        obsPos_matrix(:,i) = obs;
    end   
    cl=0;
    iters=0;
% loops through for Newton-Raphson until the norm increment is less than 10^-6
    while(iters==0)
    % define the matrix that will periodically sum the inverse of the 
    % qtilda matrices
        s_sva = length(SVA);
        inv_qtilda_all = zeros(s_sva(1)-1,s_sva(1)-1);
    % define the matrix that will periodically
    % sum the inverse times l
        inv_qtilda_lk_all = zeros(s_sva(1)-1,1);
    % define a matrix which will store all the A_k matrices
        A_k_all = [];
    % define a matrix which will store all the l_k matrices
        l_k_all = [];
    % initiates the sample count for both pieces of receiver data
        sA=0;
        sB=0;
        count=0;
        for i = 1:n
        % increments the sample
            sA = sA+1;
            sB = sB+1;
        % determine pseudo-range measurements 'pseudoR' and beat carrier
        % phase 'phi' for both receivers
            pseudoRA = pseudoA(sA,3:2:end)';
            pseudoRB = pseudoB(sB,3:2:end)';
            phiA_temp = phiA(sA,3:2:end)';
            phiB_temp = phiB(sB,3:2:end)';
        % check that have psuedo-range data for all satellites for the 
        % sample and then determine the gpsTime for each receiver
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
        % calculate qtilda, A_k, and l matrix at each time sample
        % based on the current observation guesses and environmental
        % information from that day and time
            [qtilda_inv A_k] = calculateQtilda(pseudoRA,pseudoRB,...
                    ephemA, phiA_temp,phiB_temp,gpsTimeA,...
                    obsPos_matrix(:,i),ionParam,1,elevmask,1022,9,0.86,1);
            [l_k range_known] = l_kmatrix(pseudoRA, pseudoRB,...
                    ephemA, phiA_temp, phiB_temp,gpsTimeA,gpsTimeB,...
                    obsPos_matrix(:,i),ionParam,1,elevmask,...
                    1022,9,0.86,1,pos_known);                   
        % store the current A_k matrix in a larger A_k all matrix to be 
        % used later on
            A_k_all(:,((3*i-2):3*i)) = A_k;
        % store the current l_k matrix in a larger l_k all matrix to be 
        % used later on            
            l_k_all(:,i) = l_k;
        % calculate the inverse of qtildas summed for the calculation of
        % the real-time biases
            inv_qtilda_all = inv_qtilda_all + qtilda_inv;
        % calculate the inverse of qtilda time*l's summed for the
        % calculation of the real-time biases
            qtilda_inv_l = qtilda_inv*l_k;
            inv_qtilda_lk_all = inv_qtilda_lk_all + qtilda_inv_l;
        % increment the count
            count = count + 1;
        end
    % dowload the constant file
        constant;
    % solve for the unknown biases prior to finding the exact integer 
    % solution
        biases = (1/lambdaL1)*inv(inv_qtilda_all)*inv_qtilda_lk_all; 
    % loop through and grab the correct diagonal elements for
    % sigma_graddeltabeta
        sigma_graddeltabeta_r = zeros(s_sva-1,1);
        inv_inv_qtilda_all = inv(inv_qtilda_all);
        for i = 1:s_sva-1
            sigma_graddeltabeta_r(i) = sqrt(inv_inv_qtilda_all(i,i));
        end
        jlsq_total_real = 0;
        jlsq_total_round = 0;
    % calculate aeta
        for j = 1:n
            A_temp = A_k_all(:,((3*i-2):3*i));
            l_temp = l_k_all(:,i);
        % define q_dd
            q_dd = 2.*ones(length(A_temp)) + 2.*eye(length(A_temp));
            invQdd = inv(q_dd);
        % calculate q_tilda inverse for the particular sample
            qtilda_inv_i = invQdd-(invQdd*A_temp*inv(A_temp'*invQdd*A_temp)*A_temp'*invQdd);
        % calculate jlsq
            jlsq_real = (l_temp - lambdaL1*biases)'*qtilda_inv_i*(l_temp - lambdaL1*biases);
            jlsq_round = (l_temp - lambdaL1*round(biases))'*qtilda_inv_i*(l_temp - lambdaL1*round(biases));
            jlsq_total_real = jlsq_total_real + jlsq_real;
            jlsq_total_round = jlsq_total_round + jlsq_round;
        end
        jlsq_total_real = .5*jlsq_total_real;
        jlsq_total_round = .5*jlsq_total_round;
        aeta = sqrt(.5*(jlsq_total_round-jlsq_total_real));
    % calculate the real value integer biases
        graddeltabeta_inrect = hyperrectangle(biases,...
                                               sigma_graddeltabeta_r,...
                                                aeta);
    % define a temp variable to store the lowest
    % jlsq_total and biases
        jlsq_low = 0;
        bias_low = 0;
        [x,y] = size(graddeltabeta_inrect);
    % define the outer loop for checking for best biases
    for i = 1:y
    % grab the column vector containing a possible case of the biases integer
        bias_temp = graddeltabeta_inrect(:,i);
    % define jlsq total
        jlsq_total = 0;
    % calculate jLsq(deldelBeta) in a different for loop
        for j = 1:n
            A_temp = A_k_all(:,((3*i-2):3*i));
            l_temp = l_k_all(:,i);
        %define q_dd
            q_dd = 2.*ones(length(A_temp))+2.*eye(length(A_temp));
            invQdd = inv(q_dd);
            qtilda_inv_i=invQdd-(invQdd*A_temp*inv(A_temp'*invQdd*A_temp)*A_temp'*invQdd);
        % calculate q_tilda inverse for the particular sample
        % calculate jlsq
            jlsq = (l_temp-lambdaL1*bias_temp)'*qtilda_inv_i*(l_temp-lambdaL1*bias_temp);
            jlsq_total = jlsq_total+jlsq;   
        end
        jlsq_total = .5*jlsq_total;
        if (jlsq_low == 0 || jlsq_total < jlsq_low)
            jlsq_low = jlsq_total;
            bias_low = bias_temp;   
        end    
    end                                          
    % define the biases to be the lowest cost bias
        biases = bias_low;
    % define the deltax matrix that will store all the newton raphson 
    % increments
        delta_matrix = [];
    % calculate del_x for each sample time
        for i = 1:n
        % calculate the perturbations to all the positions
            del_x=solveDDGPS(biases,A_k_all(:,((3*i-2):3*i)), l_k_all(:,i));
        % make sure del_x is a column vector
            del_x = del_x(:);
        % store the del_x value at particular sample time into matrix
            delta_matrix(:,i) = del_x;
        end
    % update the observation guess based on the newton raphson
        obsPos_matrix = obsPos_matrix + delta_matrix;
    % increment c1
        cl = cl + 1;
    %if the norm value of del_x is less then 1 stop the loop
    if(norm(del_x)<10e-6)
    iters=1;
    end
        
    end
    
    %output the surveyed and double difference gps solution of the unkown
    %reciever location
    latlong_posknown=latlong(pos_unknown_rcvr);
    Mean_xpos=mean(obsPos_matrix(1,:));
    Mean_ypos=mean(obsPos_matrix(2,:));
    Mean_zpos=mean(obsPos_matrix(3,:));
    ECEF_mean=[Mean_xpos Mean_ypos Mean_zpos];
    latlong_ecef_mean=latlong(ECEF_mean);
    fprintf('The surveyed position of the unknown reciever is: %f %f %f\n',latlong_posknown(1), latlong_posknown(2), latlong_posknown(3));
    fprintf('The average DDGPS solution of the unknown reciever is: %f %f %f\n', latlong_ecef_mean(1), latlong_ecef_mean(2), latlong_ecef_mean(3));
    %print out the normal error between the surveyed solution and
    %and the double difference gps solution
    fprintf('The norm error between the surveyed solution and the DDGPS solution is: %f\n',norm(pos_unknown_rcvr-ECEF_mean));
    %rotate the surveyed location the the double difference gps solution
    %into verticle east north by multiplying by the rotation matrices
    %R_y(-latitude)*R_z(longitude)
    degrad=pi/180;
    latOBS=latlong_posknown(1)*degrad;
    longOBS=latlong_posknown(2)*degrad;
    obsPos_ECEFx=obsPos_matrix(1,:);
    obsPos_ECEFy=obsPos_matrix(2,:);
    obsPos_ECEFz=obsPos_matrix(3,:);
    ECEFdiffx=obsPos_ECEFx-pos_unknown_rcvr(1);
    ECEFdiffy=obsPos_ECEFy-pos_unknown_rcvr(2);
    ECEFdiffZ=obsPos_ECEFz-pos_unknown_rcvr(3);
    [V_diff E_diff N_diff]=ven(latOBS,longOBS,ECEFdiffx, ECEFdiffy, ECEFdiffZ);
    %plot the East North double difference solution vs surveyed location
    %solution
    figure(1)
    plot(E_diff, N_diff, 0, 0,'g--','LineWidth',4);
    axis([-1 1 -1 1]);
    xlabel('East Antenna Seperation Error Relative to Surveyed Seperation (m)')
    ylabel('North Antenna Seperation Error Relative to Surveyed Seperation (m)')
    %plot the Verticle double difference solution vs surveyed location
    %solution vs. time
    figure(2)
    plot(pseudoA(1:n,1), V_diff);
    xlabel('Time (sec)');
    ylabel('Vertical Antenna Seperation Error Relative to Surveyed Seperation (m)');
    
    
    
    
    
    
