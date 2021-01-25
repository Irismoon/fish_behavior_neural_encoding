[sessionID,fishID] = getfish();
%first, move centerline data to local disk
% moveDataFromServer(sessionID,fishID);
%2,get angledata from centerline
% curvature_data_extract_whole_body_fish(sessionID,fishID);
%3.detect bouts
for i=7:length(sessionID)
    try
        detectBout(sessionID{i},fishID{i});
    catch ME
        disp([sessionID{i} ': ' ME.message]);
    end
    %4.select good bouts
    try
        select_bouts(sessionID{i},fishID{i});
    catch ME
        disp([sessionID{i} ': ' ME.message]);
    end
end
%%
nBout_type = zeros(length(sessionID),2);
nConv = zeros(length(sessionID),1);
for i=1:length(sessionID)
    load(fullfile(getpath('behavior',sessionID{i},fishID{i}),'tail_swing'),'good_bout_idx');
    nBout_type(i,1) = nBout_type(i,1) + nnz(good_bout_idx(:,3)==1);
    nBout_type(i,2) = nBout_type(i,2) + nnz(good_bout_idx(:,3)==2);
    nConv(i) = nnz(good_bout_idx(:,4));
    if size(good_bout_idx,1)>0
        fprintf([sessionID{i} '-' fishID{i} ': ']);
        good_bout_idx(randperm(size(good_bout_idx,1),1),:)
    end
end
figure,
iosr.statistics.boxPlot(nBout_type,'showScatter',true,'scatterMarker','o');
set(gca,'XTickLabel',{'turn','forward'});
ylabel('number of bouts');
