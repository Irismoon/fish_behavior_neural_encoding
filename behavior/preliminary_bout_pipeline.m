sessionID = '201105';fishID = '2';
%first, move centerline data to local disk
% moveDataFromServer(sessionID,fishID);
%2,get angledata from centerline 
% curvature_data_extract_whole_body_fish(sessionID,fishID);
%3.detect bouts
detectBout(sessionID,fishID);
%4.select good bouts
select_bouts(sessionID,fishID);