%% Power Law
Base =15;
mat = sparse(2^Base,2^Base);
powOut = repmat('a',Base+1,Base);
powIn = repmat('a',Base+1,Base);
edgesAdded = cell(1,Base+1);

for i = 0:Base
    % Set non-2 values different ('0' and '1') for seperate in-hub and out-hub
    % If non-2 values are same, then the in-hub and out-hub overlap
    powOut(i+1,:) = strcat(repmat('0',1,i),repmat('2',1,Base-i));
    powIn(i+1,:) = strcat(repmat('2',1,i),repmat('1',1,Base-i));
end

for i = 1:size(powOut,1)
    edgesAdded{i} = zeros(2,2^(i-1)*2^(Base+1-i));
    ops2 = Expand([powOut(i,:),powIn(i,:)]);
    [mat,edgesAdded{i}] = ConnectBarcodesSameMat(ops2(:,1:Base),ops2(:,Base+1:2*Base),mat,edgesAdded{i});
end
outDeg = full(sum(mat~=0,2));
inDeg = full(sum(mat~=0,1));
figure;
subplot(1,3,1)
spy(mat)
title('Adjacency Matrix')

% Outdegree Distribution
subplot(1,3,2)
edges = logspace(log10(min(outDeg)),log10(max(outDeg)),Base+2);
[a,b] = histcounts(outDeg,edges,'Normalization','Probability');
scatter(edges(1:Base+1),a./diff(edges));
hold on
plotx = edges(1:Base + 1);
ploty = a./diff(edges);
p = polyfit(log(plotx(1:Base)),log(ploty(1:Base)),1);
m = p(1);
b = exp(p(2));
ezplot(@(plotx) b*plotx.^m,[plotx(1) plotx(end-1)])
xlim([min(edges),max(edges)])
ylim([min(a./diff(edges)),max(a./diff(edges))])
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
xlabel('k_{out}','FontSize',14)
ylabel('p(k_{out})','FontSize',14)
title('Out Degree Distribution')
pbaspect([1 1 1])

% Indegree Distribution
subplot(1,3,3)
edges = logspace(log10(min(inDeg)),log10(max(inDeg)),Base+2);
[a,b] = histcounts(inDeg,edges,'Normalization','Probability');
scatter(edges(1:Base+1),a./diff(edges));
hold on
plotx = edges(1:Base + 1);
ploty = a./diff(edges);
p = polyfit(log(plotx(1:Base)),log(ploty(1:Base)),1);
m = p(1);
b = exp(p(2));
ezplot(@(plotx) b*plotx.^m,[plotx(1) plotx(end-1)])
xlim([min(edges),max(edges)])
ylim([min(a./diff(edges)),max(a./diff(edges))])
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
xlabel('k_{in}','FontSize',14)
ylabel('p(k_{in})','FontSize',14)
set(0,'DefaultFigureWindowStyle','docked')
title('In Degree Distribution')
pbaspect([1 1 1])

%% SF with 0.5 and 2 degree exponents
Base = 12;
mat = sparse(2^Base,2^Base);
powOut = repmat('a',Base/2,Base);
powIn = repmat('a',Base/2,Base);
edgesAdded = cell(1,Base+1);

for i = 0:Base/2
    % Set non-2 values different ('0' and '1') for seperate in-hub and out-hub
    % If non-2 values are same, then the in-hub and out-hub overlap
    powOut(i+1,:) = strcat(repmat('0',1,i),repmat('2',1,Base-i));
    powIn(i+1,:) = strcat(repmat('2',1,2*i),repmat('0',1,Base-2*i));
end

for i = 1:size(powOut,1)
    edgesAdded{i} = zeros(2,2^(i-1)*2^(Base+1-i));
    ops2 = Expand([powOut(i,:),powIn(i,:)]);
    [mat,edgesAdded{i}] = ConnectBarcodesSameMat(ops2(:,1:Base),ops2(:,Base+1:2*Base),mat,edgesAdded{i});
end
outDeg = full(sum(mat~=0,2));
inDeg = full(sum(mat~=0,1));
figure;
subplot(1,3,1)
spy(mat)
title('Adjacency Matrix')

% Outdegree Distribution
subplot(1,3,2)
r = unique(outDeg);
q = histc(outDeg, unique(outDeg))/2^Base;
scatter(r,q);
plotx = r;
hold on
p = polyfit(log10(r(1:length(r)-1)),log10(q(1:length(q)-1)),1);
m = p(1);
m
b = 10^(p(2));
ezplot(@(plotx) b*plotx.^m,[plotx(1) plotx(end)])
xlim([min(r),max(r)])
ylim([min(q),max(q)])
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
xlabel('k_{out}','FontSize',14)
ylabel('p(k_{out})','FontSize',14)
title('Out Degree Distribution')
pbaspect([1 1 1])

% Indegree Distribution
subplot(1,3,3)
r = unique(inDeg);
q = histc(inDeg, unique(inDeg))/2^Base;
scatter(r,q);
plotx = r;
hold on
p = polyfit(log10(r(1:length(r)-1)),log10(q(1:length(q)-1)),1);
m = p(1);
m
b = 10^(p(2));
b
ezplot(@(plotx) b*plotx.^m,[plotx(1) plotx(end)])
xlim([min(r),max(r)])
ylim([min(q),max(q)])
set(gca, 'xScale', 'log')
set(gca, 'yScale', 'log')
xlabel('k_{in}','FontSize',14)
ylabel('p(k_{in})','FontSize',14)
set(0,'DefaultFigureWindowStyle','docked')
title('In Degree Distribution')
pbaspect([1 1 1])

%% Binary Tree
Depth = 4;
b = Depth;
mat = sparse(2^b,2^b);
numRuleParticipate = zeros(1,2^Depth);
ToConnect = [zeros(1,Depth-1),1];

while ~isempty(ToConnect)
    CurrConnect = ToConnect(1,:);
    Rule = [CurrConnect,CurrConnect(2:Depth),2];

        ops2 = Expand(sprintf('%d', Rule));
        [mat,~] = ConnectBarcodesSameMat(ops2(:,1:Depth),ops2(:,Depth+1:2*Depth),mat,numRuleParticipate);

    if CurrConnect(2) ~= 1
        ToConnect = [ToConnect;CurrConnect(2:Depth),1;CurrConnect(2:Depth),0];
    end
    ToConnect(1,:) = [];
end
%plot(graph(symmetrize(mat(2:2^Depth,2:2^Depth))),'NodeLabel',[],'Layout','layered')
h = plot(graph(symmetrize(mat(2:2^Depth,2:2^Depth))),'Layout','layered');
labelnode(h,1:(2^b-1),cellstr(dec2bin(1:(2^b-1))))
pbaspect([1 1 1])

%% Hierarchical Network
figure 
Depth = 4;
b = (Depth+1)*3;
mat = sparse(2^b,2^b);
numRuleParticipate = zeros(1,2^b);

NonZero = [0 0 0; 1 0 0; 1 0 1; 1 1 0; 1 1 1];
Zero = [0 0 0];
Multi = [1 2 2];

dest = horzcat(reshape(NonZero(permn([1,2,3,4,5],Depth).',:).',3*Depth,5^Depth).',repmat(Multi,5^Depth,b/3-Depth));
source = horzcat(reshape(NonZero(permn([1,2,3,4,5],Depth).',:).',3*Depth,5^Depth).',repmat(Multi,5^Depth,b/3-Depth));
TotalMat = horzcat(source,dest);
for i = 1:size(TotalMat,1)
    ops2 = Expand(sprintf('%d', TotalMat(i,:)));
    [mat,~] = ConnectBarcodesSameMat(ops2(:,1:b),ops2(:,b+1:2*b),mat,numRuleParticipate);
end

while Depth > 0
    dest = horzcat(reshape(NonZero(permn([1,2,3,4,5],Depth).',:).',3*Depth,5^Depth).',repmat(Zero,5^Depth,b/3-Depth));
    source = horzcat(reshape(NonZero(permn([1,2,3,4,5],Depth).',:).',3*Depth,5^Depth).',repmat(Multi,5^Depth,b/3-Depth));
    TotalMat = horzcat(source,dest);
    for i = 1:size(TotalMat,1)
        ops2 = Expand(sprintf('%d', TotalMat(i,:)));
        [mat,~] = ConnectBarcodesSameMat(ops2(:,1:b),ops2(:,b+1:2*b),mat,numRuleParticipate);
    end
    Depth = Depth -1;
end

mat = symmetrize(mat);
mat = mat - diag(diag(mat));
[~,a] = giantComponent(mat);
G = graph(mat(a,a));


h = plot(G,'Layout','force','NodeLabel','','MarkerSize',3,'EdgeColor',[0.6,0.6,0.6],'NodeColor',[0.2,0.4,1.0]);
pbaspect([1 1 1])