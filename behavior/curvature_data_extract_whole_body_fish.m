function curvature_data_extract_whole_body_fish(sessionID,fishID)
load(fullfile(getpath('behavior',sessionID,fishID),'tail_swing'),'centerline');
    numframes=length(centerline);
    numcurvpts = 60;
    finallyPointsNum=20;
    proximity = 50;
    % chazhi to normalize centerline data num to 60
    nozerolength=zeros(numframes,1);
    for j=1:numframes
    centerline1=squeeze(centerline(j,:,:));
    a=all(centerline1==0,2);
    nozerolength(j)=140-sum(a);
    end
    wrong_tail=zeros(numframes,1);
    wrong_tail(find(nozerolength<55))=1;
    % chazhi
    new_centerline=zeros(numframes,numcurvpts,2);
    for j=1:numframes
        if(nozerolength(j)~=0)
            x=1:nozerolength(j);
            temp1=centerline(j,:,1);
            temp2=centerline(j,:,2);
            temp2(find(temp2==0))=[];
            temp1(find(temp1==0))=[];
            v1=temp1;
            v2=temp2;
            vq=1:nozerolength(j)/(numcurvpts+1):nozerolength(j)+1;
            vq1=interp1(x,v1,vq);
            vq2=interp1(x,v2,vq);
            new_centerline(j,:,1)=vq1(1:numcurvpts);
            new_centerline(j,:,2)=vq2(1:numcurvpts);   
        end
    end

    new_centerline_20=new_centerline(:,1:3:60,:);
    curvdata=zeros(numframes,finallyPointsNum);
    angledata = zeros(numframes,finallyPointsNum+1);
    Tail_position=squeeze(new_centerline_20(:,finallyPointsNum,:));

    fish_length=0;  %body length in terms of pixels

    t1=0;
    j1=0; j2=0;
    framerate=50;
    spline_p=0.005;

    for j=1:numframes     
        if(nozerolength(j)~=0)
            centerlinetmp=squeeze(new_centerline_20(j,:,:))';
            if(isnan(centerlinetmp(1,end-2)))
                centerlinetmp(:,end-2)=centerlinetmp(:,end-3);
            end
            if(isnan(centerlinetmp(1,end-1)))
                centerlinetmp(:,end-1)=centerlinetmp(:,end-2);
            end
            if(isnan(centerlinetmp(1,end)))
                centerlinetmp(:,end)=centerlinetmp(:,end-1);
            end
            new_centerline_20(j,:,:) = centerlinetmp';
        %     figure (1);
        %     plot(centerline(1,:),centerline(2,:),'k-');
        %     hold on; plot(Head_position(j,1),Head_position(j,2),'ro');
        %     hold on; plot(Tail_position(j,1),Tail_position(j,2),'bo');
        %     axis off; axis equal; hold on;
            df = diff(centerlinetmp,1,2); 
            t = cumsum([0, sqrt([1 1]*(df.*df))]);%here [0,[1:100]] adds one column by the head, thus the matrix becomes [0:101] 
            fish_length=fish_length+t(end);
            cv = csaps(t,centerlinetmp,spline_p);  

        %     figure(1);
        %     fnplt(cv, '-g'); hold off;   
        %     title(strcat(strcat(num2str(j), '/'), num2str(numframes)), 'Interpreter', 'None');
            cv2 =  fnval(cv, t)';
            df2 = diff(cv2,1,1); df2p = df2';

            splen = cumsum([0, sqrt([1 1]*(df2p.*df2p))]);
            cv2i = interp1(splen+.00001*[0:length(splen)-1],cv2, [0:(splen(end)-1)/(finallyPointsNum+1):(splen(end)-1)]);

            df2 = diff(cv2i,1,1);
            atdf2 =  unwrap(atan2(-df2(:,2), df2(:,1)));
            angledata(j,:) = atdf2';

            curv = unwrap(diff(atdf2,1)); 
            curvdata(j,:) = curv';	
        end
    end


    sum_curv=sum(curvdata(:,11:20),2);
    sum_curv(find(wrong_tail==1))=0;
    window=1000;
    sum_curv=remove_average(window,sum_curv);


    left_tail_swing=zeros(length(sum_curv),1);
    left_tail_swing(find(sum_curv>0.15))=1;
    right_tail_swing=zeros(length(sum_curv),1);
    right_tail_swing(find(sum_curv<-0.15))=-1;
    left_tail_swing(find(wrong_tail==1))=0;
    right_tail_swing(find(wrong_tail==1))=0;
    left_tailswingNum=sum(left_tail_swing);
    right_tailswingNum=-sum(right_tail_swing);

    save(fullfile(getpath('behavior',sessionID,fishID),'tail_swing.mat'),'angledata','new_centerline_20','-append');
end