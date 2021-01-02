clear all
sessid = '20201219-616125';
% sessid = '20201107-613671';
cd(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
% path = ['H:\Research\Scott\raw data\' sessid];
path = ['H:\Research\Scott\processed data\headplate\mouse'];
tseries = '1891';
% tseries = '53';
znum = '511';
% znum = '2';
rows = 'S27:U27'; % rows in Session Info to write data
% assumes Session Info is in the Dropbox\Bruker directory
% AnalyzeSession_tseries(sessid,path,tseries,znum,rows);
% cd(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
load([path '\' sessid '\RegData' sessid '-' tseries '-' znum '.mat'])
[z_reg]=Zreg_sjk1(sessid,tseries,znum);
% cd('C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');
cd('C:\Users\sjk\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');

AA = [rmsX, rmsY, z_reg];
xlswrite('Session Info.xlsx',AA,'Sheet1',rows);

%%
% z_reg=Zreg_sjk(sessid,tseries,znum);
cd('C:\Users\sjk\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');
rows = 'U28:W28'; 
AA = [rmsX, rmsY, z_reg];
xlswrite('Session Info.xlsx',AA,'Sheet1',rows);
%% to access processed data
% load(['RegData' sessid '-' tseries '.mat'])
% ppath = ['H:\Research\Scott\raw data\'];
% cd(ppath)
% 20201213-188
% Wheelspeed = readtable('20201213-616125-wheel.xlsx','Range','F1:F950');
% Time = readtable('20201213-616125-wheel.xlsx','Range','K1:K950'); 
% Time=Time-Time(1);

% Wheel = readtable('20201213-616125-wheel.xlsx','Range','F1:D950'); 
% 20201213-189
% Wheel = readtable('20201213-616125-wheel.xls','Range','D1:F871'); 
% Wheeltime=table2array(Wheel(:,1));

% Wheelspeed= table2array(Wheel(:,3));
% Wheelspeed= table2array(Wheelspeed));
figure12=figure(12);
hold on
sub(1)=subplot(2,1,1);
hold on
smx = smoothdata(abs(x_values),'gaussian',25);
plot(1:numel(smx),smx,'r','LineWidth',3);
smy = smoothdata(abs(y_values),'gaussian',25);
plot(1:numel(smy),smy,'b','LineWidth',3);
smz = smoothdata(abs(z_values),'gaussian',25);
plot(1:numel(smz),smz,'k','LineWidth',3);

xpoints = plot(1:numel(x_values),abs(x_values),'r','LineWidth',3);
xpoints.Color(4)=0.25;
ypoints = plot(1:numel(y_values),abs(y_values),'b','LineWidth',3);
ypoints.Color(4)=0.25;
zpoints=plot(1:numel(z_values),abs(z_values),'k','LineWidth',3);
zpoints.Color(4)=0.25;
% plot(1:numel(x_values),abs(x_values),'LineWidth',3);
% scatterx=scatter(1:numel(x_values),abs(x_values),30,'MarkerFaceColor','b','MarkerEdgeColor','b');
% alpha(scatterx,.2);
xtk = get(gca, 'XTick'); 
wheelinc = length(Wheelspeed)/length(xtk);
timelocs=1:wheelinc:length(Wheelspeed);
timebar = ones(length(xtk));
timelocs=round(timelocs);
    timebar = Time(timelocs);
xticklabels({timebar})
xlabel('Seconds')


ylim([0 5.5])
ylabel('Position in micrometers')
% xticklabels({'0','12','24','36','48','60','72','84','96','108','120'})
% xticklabels({'0','20','40','60','80','','82'})
xlabel('Seconds')
title('in vivo Position Measurements in Micrometers (?m) Gaussian 25');

legendPosvert = strcat('Vertical, rms = ' , num2str(rmsY),'?m');
legendPoshoriz = strcat('Horizontal, rms = ',num2str(rmsX),'?m');
legendz=strcat('Z-plane, rms = ',num2str(rmsZ),'?m');
legend(legendPosvert,legendPoshoriz,legendz);
Trim = 1:815;

sub(2)=subplot(2,1,2); hold on;
plot(Time(1:815,1),Wheelspeed(1:815,1),'LineWidth',3);
% wheelinc = length(Wheelspeed)/10;
% timelocs=1:wheelinc:length(Wheelspeed);
% timebar = ones(10);
% timelocs=round(timelocs);
%     timebar = Time(timelocs);
xtk = get(gca, 'XTick'); 
wheelinc = length(Wheelspeed(Trim))/length(xtk);
timelocs=1:wheelinc:length(Wheelspeed(Trim));
timebar = ones(length(xtk));
timelocs=round(timelocs);
    timebar = Time(timelocs);
xticklabels({timebar})
xlabel('Seconds')
ylabel('Mouse Walking Velocity (cm/s)');
% linkaxes(sub,'x');
%%
% makefigs
% sum_plot

%%
sessid = '20201107-613671';
path = ['D:\Research\Scott\raw data\' sessid];
tseriesrange = ['']; % t series you want to analyze
znum = '2'; % which z-stack to use as the reference; i usually do the single
            % first block) to first assess which z series is the most
            % precise
rowrange = ['' ''];
for ii =tseriesrange(1):tseriesrange(end)
    tseries = tseriesrange(ii);
    
    for jj =rowrange(1):rowrange(end)
    rows = ['S' rowrange(jj) ':U' rowrange(jj)]; % rows in Session Info to write data
    
    AnalyzeSession_tseries(sessid,path,tseries,znum,rows);
        % assumes Session Info is in the Dropbox\Bruker directory
        % this will autofil the excel sheet Session Info
    end
end
%%
cd('C:\Users\sjk\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Figures\experimental');
%    
   z_values = thisshiftgauss-thisshiftgauss(2);
       
   Z=z_values-mean(z_values);
   arithmetic_mean_of_squares=sum(Z.^2);
   rmsZ=sqrt(arithmetic_mean_of_squares/length(z_values));
%displacement
for qq=1:length(x_values)-1
    deltax2(qq) = x_values(qq+1)-x_values(qq);
end
for ww=1:length(y_values)-1
    deltay2(ww)=y_values(ww+1)-y_values(ww);
end
for ee=1:length(z_values)-1
    deltaz2(ee)=z_values(ee+1)-z_values(ee);
end
        X=deltax2-mean(deltax2);
    Y=deltay2-mean(deltay2);
    Z=deltaz2-mean(deltaz2);
    arithmetic_mean_of_squares=sum(X.^2);
    delrmsX=sqrt(arithmetic_mean_of_squares/length(deltax2));
    arithmetic_mean_of_squares=sum(Y.^2);
    delrmsY=sqrt(arithmetic_mean_of_squares/length(deltay2));
    
    arithmetic_mean_of_squares=sum(Z.^2);
    delrmsZ=sqrt(arithmetic_mean_of_squares/length(deltaz2));


    
figure8 = figure(8);
    
rx = 1:numel(x_values);
ry = [1:length(y_values)];
rz=[1:length(z_values)];
hold on
scatterx=scatter(rx,abs(x_values),65,'MarkerFaceColor','r','MarkerEdgeColor','k');
alpha(scatterx,.3);
scattery=scatter(ry,abs(y_values),65,'MarkerFaceColor','b','MarkerEdgeColor','k');
alpha(scattery,.3);
scatterz=scatter(rz,abs(z_values),65,'MarkerFaceColor','k','MarkerEdgeColor','k');
alpha(scatterz,.3);


ylim([0 5])
ylabel('Position in micrometers')
xticklabels({'0','12','24','36','48','60','72','84','96','108','120'})
xlabel('Seconds')
title('in vivo Position Measurements (?m)');
legendPosvert = strcat('Vertical, rms = ' , num2str(rmsY),'?m');
legendPoshoriz = strcat('Horizontal, rms = ',num2str(rmsX),'?m');
legendz=strcat('Z-plane, rms = ',num2str(rmsZ),'?m');
legend(legendPosvert,legendPoshoriz,legendz);

saveas(figure8,'20201213-616125-188-50-position-xyz.eps')

%

figure9=figure(9);
hold on
scatterx=scatter(1:numel(deltax2),abs(deltax2),65,'MarkerFaceColor','r','MarkerEdgeColor','k');
alpha(scatterx,.3);
scattery=scatter(1:numel(deltay2),abs(deltay2),65,'MarkerFaceColor','b','MarkerEdgeColor','k');
alpha(scattery,.3);
scatterz=scatter(1:numel(deltaz2),abs(deltaz2),65,'MarkerFaceColor','k','MarkerEdgeColor','k');
alpha(scatterz,.3);


ylim([0 5.5])
ylabel('Displacement in micrometers')
xticklabels({'0','12','24','36','48','60','72','84','96','108','120'})
xlabel('Seconds')
title('in vivo Displacement Measurements (?m)');

legendPosvert = strcat('Vertical, rms = ' , num2str(delrmsY),'?m');
legendPoshoriz = strcat('Horizontal, rms = ',num2str(delrmsX),'?m');
legendz=strcat('Z-plane, rms = ',num2str(delrmsZ),'?m');
legend(legendPosvert,legendPoshoriz,legendz);
saveas(figure9,'20201213-616125-188-50-displacement-xyz.eps')
%%
fig10 = figure(10);
hold on
smx = smoothdata(abs(x_values),'gaussian',25);
plot(1:numel(smx),smx,'r','LineWidth',3);
smy = smoothdata(abs(y_values),'gaussian',25);
plot(1:numel(smy),smy,'b','LineWidth',3);
smz = smoothdata(abs(z_values),'gaussian',25);
plot(1:numel(smz),smz,'k','LineWidth',3);

xpoints = plot(1:numel(x_values),abs(x_values),'r','LineWidth',3);
xpoints.Color(4)=0.25;
ypoints = plot(1:numel(y_values),abs(y_values),'b','LineWidth',3);
ypoints.Color(4)=0.25;
zpoints=plot(1:numel(z_values),abs(z_values),'k','LineWidth',3);
zpoints.Color(4)=0.25;
% plot(1:numel(x_values),abs(x_values),'LineWidth',3);
% scatterx=scatter(1:numel(x_values),abs(x_values),30,'MarkerFaceColor','b','MarkerEdgeColor','b');
% alpha(scatterx,.2);


ylim([0 5.5])
ylabel('Position in micrometers')
xticklabels({'0','12','24','36','48','60','72','84','96','108','120'})
xlabel('Seconds')
title('in vivo Position Measurements in Micrometers (?m) Gaussian 25');

legendPosvert = strcat('Vertical, rms = ' , num2str(rmsY),'?m');
legendPoshoriz = strcat('Horizontal, rms = ',num2str(rmsX),'?m');
legendz=strcat('Z-plane, rms = ',num2str(rmsZ),'?m');
legend(legendPosvert,legendPoshoriz,legendz);
% saveas(fig10,'1213-616125-188-50-posgaus25.eps')