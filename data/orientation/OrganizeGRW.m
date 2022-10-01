% OrganizeGRW takes the data from the study "Accuracy and response times 
% in a continuous task" and organizes it into individual variables
% corresponding to response time, response deviation, manipulations /
% factor levels for each trial.

% Author: Peter D. Kvam
% Date: August 30, 2016
% Contact the author with questions at kvam.peter@gmail.com

%% Move data from individual files to vectors
fileExts = [100,101,110,120,130,140,150,170,180,190,200,210]; % names of files containing data

trialIndex = 0;


for s = fileExts
    thisFile = ['GRW',num2str(s),'.mat'];
    load(thisFile);
    nTrials = length(dataMat.trial);    

    for n = 1:nTrials
        trialIndex = trialIndex + 1;
        
        sNum(trialIndex) = s;
        target(trialIndex) = dataMat.trial(n).target;
        resp(trialIndex) = dataMat.trial(n).resp;
        dev(trialIndex) = dataMat.trial(n).dev;
        rt(trialIndex) = dataMat.trial(n).RT;
        isSpeed(trialIndex) = dataMat.trial(n).speedCond;
        isCued(trialIndex) = dataMat.trial(n).cueCond;
        cueOrient(trialIndex) = dataMat.trial(n).cueOrient;
        points(trialIndex) = dataMat.trial(n).trialPoints;
        jitter(trialIndex) = dataMat.trial(n).jitter;
    end
end

%% Clean outliers

notBadRTs = find((rt < 2.5) .* (rt > .15));

target = target(notBadRTs);
resp = resp(notBadRTs);
dev = dev(notBadRTs);
rt = rt(notBadRTs);
isSpeed = isSpeed(notBadRTs);
isCued = isCued(notBadRTs);
cueOrient = cueOrient(notBadRTs);
points = points(notBadRTs);
jitter = jitter(notBadRTs);
    
          
% Calculate absolute deviations and find list of indices for speed & accuracy conditions
absDev = abs(dev);
speedTrials = find(isSpeed);
accTrials = find(~isSpeed);



% Get cue orientation, rt, and deflections for cued condition
cueOrientCued = cueOrient(find(isCued));
rtCued = rt(find(isCued));
accCued = absDev(find(isCued)) .* 180 / pi;


cueDefs(find(isCued)) = target(find(isCued)) - cueOrientCued; % deflection of cue from actual target
cueDefs(find(cueDefs > (pi/2))) = cueDefs(find(cueDefs > (pi/2))) - pi;
cueDefs = abs(round(cueDefs .* 100 .* 180 ./ pi)./100) ;
% cueDefs = round(100*mod(cueDefs + pi, pi/2))/100; % put deflections onto [0,pi]
[uCues,u2,cueIndex] = unique(cueDefs); % u3 will correspond to the index of the cue orientation in u1