function codes = Expand(twoVec)
twoLocs = find(twoVec=='2');
numTwos = length(twoLocs);
ops = dec2base(0:2^numTwos-1,2);
codes = repmat(twoVec,[length(ops),1]);
codes(:,twoLocs) = ops;