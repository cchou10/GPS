%a function which runs pseudorangeerrors in order
%to determine the neccesary time samples
pseudorangeerrors;
prerrhist=prerrorPRN01;
npolyfit=1;
[ Corrvec, mlagsvec, polyprerr, prerrresidhist] = ...
                xcorrpolyfit(prerrhist,npolyfit);
plot(mlagsvec,Corrvec);