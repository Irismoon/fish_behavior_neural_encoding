%%
%investigate the raw trace of region #356 which contains most 2nd PC time
%series
% load('tail_swing','bout_idx');

spk_focus = Spike_X_EstTrace(356,:);
figure,
group = zeros(numWholeTime,1);group(bout_idx_aligned_wholecourse) = 1;group(nonbout_idx)=2;
subplot(2,1,1)
p = ranksum(spk_focus(bout_idx_aligned_wholecourse),spk_focus(nonbout_idx));
boxplot(spk_focus,group);
sigstar([1 2],p);
set(gca,'XTickLabel',{'bout','nonbout'});
title(['Comparison of neural activity between bout and no bouts: p=' num2str(p)]);
%control
spk_focus_ctrl = spk_focus(randperm(length(spk_focus)));
%and rerun the code above
subplot(2,1,2),
p = ranksum(spk_focus_ctrl(bout_idx_aligned_wholecourse),spk_focus_ctrl(nonbout_idx));
boxplot(spk_focus_ctrl,group);
sigstar([1 2],p);
set(gca,'XTickLabel',{'bout','nonbout'});
title(['control: p=' num2str(p)]);
p = ranksum(spk_focus(hunting_idx_wholecourse),spk_focus(setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)));
figure,
subplot(2,1,1),
group = 2*ones(length(bout_idx_aligned_wholecourse),1);
group(ismember(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)) = 1;
boxplot(spk_focus(bout_idx_aligned_wholecourse),group);
hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'hunting bouts','nonhunting bouts'});
subplot(2,1,2),
group = 2*ones(length(bout_idx_aligned_wholecourse),1);
group(ismember(bout_idx_aligned_wholecourse,hunting_idx_wholecourse)) = 1;
boxplot(spk_focus(bout_idx_aligned_wholecourse),group(randperm(length(bout_idx_aligned_wholecourse))));
hold on;
idxytmp=randperm(length(bout_idx_aligned_wholecourse),length(hunting_idx_wholecourse));
p = ranksum(spk_focus(bout_idx_aligned_wholecourse(idxytmp)),spk_focus(bout_idx_aligned_wholecourse(setdiff(1:length(bout_idx_aligned_wholecourse),idxytmp))));
sigstar([1 2],p);
set(gca,'XTickLabel',{'hunting bouts','nonhunting bouts'});
title('control');
sgtitle('comparison between hunting and nonhuting bouts');
set(gca,'XTickLabel',{'hunting','nonhunting'});
%compare nonhunting bouts and nonbouts
figure,
x1 = spk_focus(setdiff(1:numWholeTime,bout_idx_aligned_wholecourse))';
x2 = spk_focus(setdiff(bout_idx_aligned_wholecourse,hunting_idx_wholecourse))';
x = cat(1,x2,x1);
group = [ones(length(x2),1);2*ones(length(x1),1)];
p = ranksum(x1,x2);
subplot(2,1,1),
boxplot(x,group);hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'nonhunting bouts','nonbouts'});
subplot(2,1,2),
randidx = randperm(length(x));
p = ranksum(x(randidx(1:length(x1))),x(randidx((length(x1)+1):end)));
boxplot(x,group(randidx));hold on;
sigstar([1 2],p);
sgtitle('comparison between nonhunting bouts and nonbouts');
%compare left and right bouts
figure,
x1 = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',find(param_head_angle_hunting>0),'un',0);
x1 = spk_focus(cat(1,x1{:}))';%left
x2 = arrayfun(@(i) (hunting_idx_start(i):hunting_idx_end(i))',find(param_head_angle_hunting<0),'un',0);
x2 = spk_focus(cat(1,x2{:}))';%right
x = cat(1,x2,x1);
group = [ones(length(x2),1);2*ones(length(x1),1)];
p = ranksum(x1,x2);
subplot(2,1,1),
boxplot(x,group);hold on;
sigstar([1 2],p);
set(gca,'XTickLabel',{'right','left'});
subplot(2,1,2),
randidx = randperm(length(x));
p = ranksum(x(randidx(1:length(x1))),x(randidx((length(x1)+1):end)));
boxplot(x,group(randidx));hold on;
sigstar([1 2],p);
sgtitle('comparison between left and right hunting bouts');
%%
%compare left and right vision
x1 = spk_focus(hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)>0))';
x2 = spk_focus(hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)<0))';
px = ranksum(x1,x2);

y1 = spk_focus(simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)>0))';
y2 = spk_focus(simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)<0))';
y = cat(1,y1,y2);
py = ranksum(y1,y2);

z1 = spk_focus(rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)>0))';
z2 = spk_focus(rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)<0))';
pz = ranksum(z1,z2);
data = {x1,x2;y1,y2;z1,z2};
h = boxplot2(data);
set(gca,'XTick',[1 2 3],'XTickLabel',{'hunting','simultaneous','resting'});
sigstar([h.lwhis(1,1).XData(1) h.lwhis(2,1).XData(1)],px);
sigstar([h.lwhis(1,2).XData(1) h.lwhis(2,2).XData(1)],py);
sigstar([h.lwhis(1,3).XData(1) h.lwhis(2,3).XData(1)],pz);
sigstar([h.lwhis(1,1).XData(1) h.lwhis(1,2).XData(1)],ranksum(x1,y1));
sigstar([h.lwhis(1,2).XData(1) h.lwhis(1,3).XData(1)],ranksum(z1,y1));
sigstar([h.lwhis(2,1).XData(1) h.lwhis(2,2).XData(1)],ranksum(x2,y2));
sigstar([h.lwhis(2,2).XData(1) h.lwhis(2,3).XData(1)],ranksum(z2,y2));
title('Left Right visual encoding');
%%
%if the left/right difference comes from the distance/angle
idx_no_param = find(param_head_dist_fluo>300);
hunting_idx_wholecourse(ismember(hunting_idx_wholecourse,idx_no_param)) = [];
simul_bout_idx_wholecourse(ismember(simul_bout_idx_wholecourse,idx_no_param)) = [];
rest_idx_wholecourse(ismember(rest_idx_wholecourse,idx_no_param)) = [];
dx = param_head_dist_fluo(hunting_idx_wholecourse);
dy = param_head_dist_fluo(simul_bout_idx_wholecourse);
dz = param_head_dist_fluo(rest_idx_wholecourse);
data = {dx,dy,dz};
figure,
subplot(1,2,1),
boxplot2(data');
pxy = ranksum(dx,dy);pyz = ranksum(dy,dz);
hold on;
sigstar([1 2],pxy);sigstar([2 3],pyz);
ylabel('param distance');
set(gca,'XTick',[1 2 3],'XTickLabel',{'hunitng','simultaneous','resting'});
subplot(1,2,2),
anglex = param_head_angle_fluo(hunting_idx_wholecourse);
angley = param_head_angle_fluo(simul_bout_idx_wholecourse);
anglez = param_head_angle_fluo(rest_idx_wholecourse);
data = {anglex,angley,anglez};
boxplot2(data');

%%
%behavior state x left/right
behav_var = param_head_dist_fluo;
idx_behav7direct = cell(3,2);
idx_behav7direct{1,1} = hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)>0);
idx_behav7direct{1,2} = hunting_idx_wholecourse(param_head_angle_fluo(hunting_idx_wholecourse)<0);
idx_behav7direct{2,1} = simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)>0);
idx_behav7direct{2,2} = simul_bout_idx_wholecourse(param_head_angle_fluo(simul_bout_idx_wholecourse)<0);
idx_behav7direct{3,1} = rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)>0);
idx_behav7direct{3,2} = rest_idx_wholecourse(param_head_angle_fluo(rest_idx_wholecourse)<0);
spk_behav7direct = cellfun(@(x) spk_focus(x)',idx_behav7direct,'un',0);
dist_behav7direct = cellfun(@(x) behav_var(x),idx_behav7direct,'un',0);
xlim = [min(cellfun(@min,dist_behav7direct),[],'all') max(cellfun(@max,dist_behav7direct),[],'all')];
ylim = [min(cellfun(@min,spk_behav7direct),[],'all') max(cellfun(@max,spk_behav7direct),[],'all')];
figure,
coeff = zeros(3,2,2);GOF = zeros(3,2);
Edge = cell(3,2);M_SPK = cell(3,2);Weight = cell(3,2);
for i=1:3
    for j=1:2
        subplot(3,2,(i-1)*2+j),
%         scatter(dist_behav7direct{i,j},spk_behav7direct{i,j},10,'filled');hold on;
        [Y,edges] = discretize(dist_behav7direct{i,j},0:1:max(dist_behav7direct{i,j}));
        m_spk = arrayfun(@(k) quantile(spk_behav7direct{i,j}(Y==k),0.75),1:(length(edges)-1));
%         m_spk(isnan(m_spk)) = 0;        
        [N,xedge] = histcounts(dist_behav7direct{i,j},edges);
        m_spk(N==0) = [];
        edges(N==0) = [];
        xedge(N==0) = [];
        N(N==0) = [];
        Edge{i,j} = edges(1:end-1);
        M_SPK{i,j} = m_spk;
        Weight{i,j} = N;
%         scatter(edges(1:end-1),m_spk,N);
        [f,gof] = fit(edges(1:end-1)',m_spk','poly1','weights',N/sum(N));
        coeff(i,j,:) = [f.p1 f.p2];% f.p3];
        GOF(i,j) = gof.rsquare;
%         yest = f.p1*edges(1:end-1)+f.p2;
%         hold on;plot(edges(1:end-1),yest);
%         text(edges(1),yest(1),num2str(gof.rsquare));
%         p_obs = f.p1;
%         p_shuffle = zeros(1000,1);
        %shuffle test
%         for iperm = 1:1000
%             randidx = randperm(length(m_spk));
%             f = fit(edges(1:end-1)',m_spk(randidx)','poly1','weights',N/sum(N));
%             p_shuffle(iperm) = f.p1;
%         end
%         h = ttest(p_shuffle)
%         title(['normal? :' num2str(h) 'p value:' num2str(mean(p_shuffle<p_obs),3)]);
%         [f,gof] = fit(edges(1:end-1)',m_spk','exp1','weights',N/sum(N));
%         yest =f.a*exp(f.b*edges(1:end-1));
%         hold on;plot(edges(1:end-1),yest);
%         text(edges(end-1),yest(end-1),num2str(gof.rsquare));
%         yyaxis right;
%         plot(xedge(1:end-1),N);
%         axis tight;
%         [N,xedge,yedge] = histcounts2(dist_behav7direct{i,j},spk_behav7direct{i,j},'BinWidth',[2 0.01],'Normalization','probability','XBinLimits',xlim,'YbinLimits',ylim);
%         N = N./sum(N,2);
%         imagesc(N);
%         set(gca,'XTick',1:length(yedge),'XTickLabel',yedge,'YTick',1:length(xedge),'YTickLabel',xedge);
    end
end
xlabel('distance');ylabel('probability');yyaxis left;ylabel('median firing rate');

%demix the distance and behavior state condition
figure,
color = rgb({'darkgreen','magenta','navy'});
hold on,
lateral = 2;
for i=1:3
   scatter(Edge{i,lateral},M_SPK{i,lateral},Weight{i,lateral},color(i,:),'filled');
%    h(i) = plot(Edge{i,lateral},coeff(i,lateral,1)*(Edge{i,lateral}.^2)+coeff(i,lateral,2)*Edge{i,lateral}+coeff(i,lateral,3),'Color',color(i,:));
    h(i) = plot(Edge{i,lateral},coeff(i,lateral,1)*Edge{i,lateral}+coeff(i,lateral,2),'Color',color(i,:));
    text(Edge{i,lateral}(1),M_SPK{i,lateral}(1),num2str(GOF(i,lateral)),'Color',color(i,:));
end
legend(h,{'hunting','simultaneous','resting'});
%%
 %select the timeponits of simultaneous state which has the same distance
 %with hunting, and compare their firing rates
 idxxtmp = find(anglex<0);
 dx_right = dx(idxxtmp);
dlim = [min(dx_right) quantile(dx_right,0.75)];
idxytmp = find(dy<dlim(2)&dy>dlim(1)&angley<0);
idxytmp = idxytmp(randperm(length(idxytmp),length(dx_right)));
dy_sel = dy(idxytmp);
spk_y_sel = spk_focus(simul_bout_idx_wholecourse(idxytmp));
[a,b] = min(abs(dy_sel-dx_right'),[],2);
thld = 0.2;
dy_sam = dy_sel(a<thld);
pair = [find(a<thld) b(a<thld)];
figure,
subplot(3,1,1),
boxplot([dy_sel dx_right],'notch','on');
p = ranksum(dx_right,dy_sel);
hold on;sigstar([1 2],p);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[dy_sel(pair(i,1)) dx_right(pair(i,2))],'--ok','MarkerSize',3),1:size(pair,1),'un',0);
ylabel('distance');
subplot(3,1,2),
spk_x = spk_focus(hunting_idx_wholecourse);
spk_x = spk_x(idxxtmp);
boxplot([spk_y_sel spk_x],'notch','on');
p = ranksum(spk_y_sel,spk_x);
hold on;
sigstar([1 2],p);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[spk_y_sel(pair(i,1)) spk_x(pair(i,2))],'--ok','MarkerSize',3),1:size(pair,1),'un',0);
set(gca,'XTickLabel',{'randomly selected simultaneous','hunting'});
ylabel('Firing Rate');
sgtitle('Compare firing rate at same distance but diff behavior');
subplot(3,1,3),
boxplot([angley(idxytmp) anglex(idxxtmp)]);
hold on;
p = ranksum(angley(idxytmp),anglex(idxxtmp));
sigstar([1 2],p);
incre = arrayfun(@(i) spk_x(pair(i,2))-spk_y_sel(pair(i,1)),1:size(pair,1));
cmap = colormap(jet(length(incre)));
[~,I] = sort(incre,'ascend');
cmap = cmap(I,:);
angley_sel = angley(idxytmp);
arrayfun(@(i) plot([1+randn*0.03 2+randn*0.03],[angley_sel(pair(i,1)) anglex(idxxtmp(pair(i,2)))],'o-','MarkerSize',3,'Color',cmap(i,:)),1:size(pair,1),'un',0);
ylabel('angle');
%%
figure,
% polarscatter(deg2rad(anglex),ones(length(anglex),1)+randn(length(anglex),1)*0.05,5,spk_focus(hunting_idx_wholecourse),'filled');colormap('jet');
% polarscatter(deg2rad(angley),ones(length(angley),1)+randn(length(angley),1)*0.05,5,spk_focus(simul_bout_idx_wholecourse),'filled');colormap('jet');
% polarscatter(deg2rad(anglez),ones(length(anglez),1)+randn(length(anglez),1)*0.05,5,spk_focus(rest_idx_wholecourse),'filled');colormap('jet');
spkx = spk_focus(hunting_idx_wholecourse);
spky = spk_focus(simul_bout_idx_wholecourse);
spkz = spk_focus(rest_idx_wholecourse);
spk3 = cat(1,spkx,spky,spkz);
angle3 = cat(1,anglex,angley,anglez);
spkx_right = spkx(anglex<0);
spky_right = spky(angley<0);
spkz_right = spkz(anglez<0);
spk3_right = cat(1,spkx_right,spky_right,spkz_right);
angle3_right = cat(1,anglex(anglex<0),angley(angley<0),anglez(anglez<0));
[N,Xedges,Yedges] = histcounts2(spk3,angle3,'BinWidth',[0.02 2]);
% imagesc(N./sum(N,1));
% imagesc(N);
spkm = arrayfun(@(i) quantile(spk3(angle3<Yedges(i+1)&angle3>Yedges(i)),.75),1:(length(Yedges)-1));
scatter(Yedges(1:end-1),spkm,10,'filled');
% xlabel('angle');ylabel('firingRate');
% set(gca,'XTick',5:5:length(Xedges),'XTickLabel',Xedges(5:5:length(Xedges)),'YTick',5:5:length(Yedges),'YTickLabel',Yedges(5:5:length(Yedges)));
%%
%comprehensively consider distance, angle at right side, and behavior
%states and see which is important for the neural activity
dx_right = dx(anglex<0);
dy_right = dy(angley<0);
dz_right = dz(anglez<0);
d3_right = cat(1,dx_right,dy_right,dz_right);
behavior_states = [repmat({'hunting'},length(dx_right),1);repmat({'simultaneous'},length(dy_right),1);repmat({'rest'},length(dz_right),1)];
[p,tb1] = anovan(spk3_right,{behavior_states,d3_right,angle3_right},'continuous',[2 3],'model','interaction','varnames',{'behavior states','distance','angle'});
%%
%consider the relative angle between prey and eye rather than that between
%the prey and head
leftanglex = left_eye_to_param_angle_fluo(hunting_idx_wholecourse);
rightanglex = right_eye_to_param_angle_fluo(hunting_idx_wholecourse);
leftangley = left_eye_to_param_angle_fluo(simul_bout_idx_wholecourse);
rightangley = right_eye_to_param_angle_fluo(simul_bout_idx_wholecourse);
leftanglez = left_eye_to_param_angle_fluo(rest_idx_wholecourse);
rightanglez = right_eye_to_param_angle_fluo(rest_idx_wholecourse);
leftanglex_right = leftanglex(anglex<0);
rightanglex_right = rightanglex(anglex<0);
leftangley_right = leftangley(angley<0);
rightangley_right = rightangley(angley<0);
leftanglez_right = leftanglez(anglez<0);
rightanglez_right = rightanglez(anglez<0);
leftangle3_right = cat(1,leftanglex_right,leftangley_right,leftanglez_right);
rightangle3_right = cat(1,rightanglex_right,rightangley_right,rightanglez_right);
[p,tb2] = anovan(spk3_right,{behavior_states,d3_right,leftangle3_right,rightangle3_right},'continuous',[2 3 4],'model','interaction','varnames',{'behavior states','distance','leftEyeAngle','rightEyeAngle'});
%%
param_speedx = param_speed_fluo(hunting_idx_wholecourse);
param_speedy = param_speed_fluo(simul_bout_idx_wholecourse);
param_speedz = param_speed_fluo(rest_idx_wholecourse);
param_speedx_right = param_speedx(anglex<0);
param_speedy_right = param_speedy(angley<0);
param_speedz_right = param_speedz(anglez<0);
param_speed3_right = cat(1,param_speedx_right,param_speedy_right,param_speedz_right);
[p,tb2] = anovan(spk3_right,{behavior_states,d3_right,leftangle3_right,rightangle3_right,param_speed3_right},...
    'continuous',[2 3 4 5],'model','interaction','varnames',...
    {'behavior states','distance','leftEyeAngle','rightEyeAngle','speed'});
%%
%the movement speed and angle relative to the right eye
% param_righteye_pos = param_righteye_dist.*[cos(right_eye_to_param_angle) sin(right_eye_to_param_angle)];
param_righteye_pos = [param_pos_all(:,2) param_pos_all(:,1)];
param_righteye_vel = diff(param_righteye_pos,1,1);
param_righteye_speed = [0;sqrt(param_righteye_vel(:,1).^2 + param_righteye_vel(:,2).^2)];
param_righteye_move_angle = [0;rad2deg(atan(param_righteye_vel(:,2)./param_righteye_vel(:,1)))];
param_righteye_speed = param_righteye_speed(align_with_fluo_low==1);
param_righteye_move_angle = param_righteye_move_angle(align_with_fluo_low==1);
param_righteye_speed_fluo = param_righteye_speed(mask);
param_righteye_move_angle_fluo = param_righteye_move_angle(mask);
param_re_speedx_right = param_righteye_speed_fluo(hunting_idx_wholecourse(anglex<0));
param_re_speedy_right = param_righteye_speed_fluo(simul_bout_idx_wholecourse(angley<0));
param_re_speedz_right = param_righteye_speed_fluo(rest_idx_wholecourse(anglez<0));
param_re_speed3_right = cat(1,param_re_speedx_right,param_re_speedy_right,param_re_speedz_right);
param_re_moveAnglex_right = param_righteye_move_angle_fluo(hunting_idx_wholecourse(anglex<0));
param_re_moveAngley_right = param_righteye_move_angle_fluo(simul_bout_idx_wholecourse(angley<0));
param_re_moveAnglez_right = param_righteye_move_angle_fluo(rest_idx_wholecourse(anglez<0));
param_re_moveAngle3_right = cat(1,param_re_moveAnglex_right,param_re_moveAngley_right,param_re_moveAnglez_right);
[p,tb1] = anovan(spk3_right,{param_re_moveAngle3_right,param_re_speed3_right},'continuous',[1 2],'model','interaction','varname',{'velocity angle','param_speed'});
%%
%decode left/right,distance,righteyeangle using time points at different behavior states
% spk3 = {spkx;spky;spkz};
spk3 = {spkx_right;spky_right;spkz_right};
d3 = {dx_right;dy_right;dz_right};
righteyeAngle3 = {rightanglex_right;rightangley_right;rightanglez_right};
% angle3 = {anglex;angley;anglez};
angle3 = {anglex(anglex<0);angley(angley<0);anglez(anglez<0)};
eval_d3 = zeros(3,1);eval_righteyeAngle = zeros(3,1);eval_angle = zeros(3,1);
for i=1:3
   [f,gof] = fit(spk3{i},d3{i},'poly1');
   eval_d3(i) = gof.rsquare;
   [f,gof] = fit(spk3{i},cos(deg2rad(righteyeAngle3{i})),'poly1');
   eval_righteyeAngle(i) = gof.rsquare;
   [f,gof] = fit(spk3{i},cos(deg2rad(angle3{i})),'poly1');
   eval_angle(i) = gof.rsquare;
end