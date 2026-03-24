% MakeColorBlurFigures.m
%
% Generates figures 2, 3 and 4 for the color blur paper.  

allData = xls2struct('all_data.xlsx');

pcThresh = 0.75;  % should match the paper's selection criterion
colList = [.25,.25,.25;1,0,0];

% Figure 2A.   experiment 3 subject-subject plot of percent correct

cx = sin(linspace(-pi,pi,101));
cy = cos(linspace(-pi,pi,101));

exp3Data = subStruct(allData,strcmp(allData.exp,'Exp3'));
subs  = unique(exp3Data.subId);
nSubs = length(subs);

pc = zeros(nSubs,2);
rt = zeros(nSubs,2);
for i=1:length(subs)
    thisSub3bw = subStruct(exp3Data,strcmp(exp3Data.subId,subs{i})' & [exp3Data.color{:}] == 0);
    pc(i,1) = 100*mean([thisSub3bw.correct{:}]);
    rt(i,1) = mean([thisSub3bw.rt{:}]);

    thisSub3color = subStruct(exp3Data, strcmp(exp3Data.subId,subs{i})' & [exp3Data.color{:}] == 1);
    pc(i,2) = 100*mean([thisSub3color.correct{:}]);
    rt(i,2) = mean([thisSub3color.rt{:}]);
end

figure(1)
clf
hold on
rad =1;
plot([50,100],[50,100],'k-')
for i=1:nSubs
    patch(pc(i,1)+rad*cx,pc(i,2)+rad*cy,'b','FaceAlpha','.5')
end

mpc = mean(pc/100);
sem = 100*sqrt(mpc.*(1-mpc)./nSubs);
m = mean(pc);

plot(m(1),m(2),'ws','MarkerSize',12,'MarkerFaceColor','k')

axis equal
set(gca,'XLim',[50,100])
set(gca,'YLim',[50,100])
set(gca,'XTick',50:10:100)
set(gca,'YTick',50:10:100)
grid
xlabel('Grayscale Percent Correct')
ylabel('Color Percent Correct')

%%
% Figure 2B, experiment 3 subject-subject plot of RT's

rtRange = [0,3.5];

figure(2)
clf
hold on
rad =.07;
plot(rtRange,rtRange,'k-')
for i=1:nSubs
    patch(rt(i,1)+rad*cx,rt(i,2)+rad*cy,'b','FaceAlpha','.5')
end

mrt = mean(rt);
sem = std(rt)/sqrt(size(rt,1));

plot(mrt(1),mrt(2),'ws','MarkerSize',12,'MarkerFaceColor','k')



axis equal
axis tight
set(gca,'XLim',rtRange)
set(gca,'YLim',rtRange)
% set(gca,'XTick',50:10:100)
% set(gca,'YTick',50:10:100)
grid
xlabel('Grayscale RT (s)')
ylabel('Color RT (s)')


%%
%  Figures 3 and 4, Percent correct across blur levels for color and
%  grayscale

expNameList = {'Exp1','Exp2'};  % original naming conventions

for e = 1:length(expNameList)

    expName = expNameList{e}; %'Exp1' or 'Exp2' (original naming conventions)

    subData = subStruct(allData,strcmp(allData.exp,expName)' & [allData.pc{:}]>pcThresh);
    subs  = unique(subData.subId);
    blurList =unique([subData.diopter{:}]);
    nSubs = length(subs);
    nBlurs = length(blurList);
    c = zeros(nBlurs,2,nSubs);
    n = zeros(nBlurs,2,nSubs);

    for i=1:nSubs
        for j=1:nBlurs
            thisSubbw = subStruct(subData,strcmp(subData.subId,subs{i})' & ...
                [subData.diopter{:}] == blurList(j) & ...
                [subData.color{:}] == 0);
            c(j,1,i) =  sum([thisSubbw.correct{:}]);
            n(j,1,i) = length([thisSubbw.correct{:}]);

            thisSubcolor = subStruct(subData,strcmp(subData.subId,subs{i})' & ...
                [subData.diopter{:}] == blurList(j) & ...
                [subData.color{:}] == 1);
            c(j,2,i) =  sum([thisSubcolor.correct{:}]);
            n(j,2,i) = length([thisSubcolor.correct{:}]);
        end
    end


    pc = mean(c./n,3);

    sem = sqrt(pc.*(1-pc)./nSubs);
    figure(e+2)
    clf
    hold on
    plot(blurList,50*ones(size(blurList)),'k:')
    for i=1:2
        errorbar(blurList,100*pc(:,i),100*sem(:,i),'Color',colList(i,:),'LineStyle','none')
        h(i) = plot(blurList,100*pc(:,i),'o-','Color',colList(i,:),'MarkerFaceColor',colList(i,:));
    end
    set(gca,'XLim',[-.25,max(blurList)+.25])

    legend(h,{'Grayscale','Color'},'Location','NorthEast')
    ylabel('Percent Correct')
    xlabel('Diopter')
    title(expName)
end

%%
% support function

function subDat = subStruct(dat,id)

subDat = struct();
fields = fieldnames(dat);
for i=1:length(fields)
    subDat.(fields{i}) = dat.(fields{i})(id);
end

end
