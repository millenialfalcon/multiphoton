clear all
sessid = '20201228-616125';
% sessid = '20201107-613671';
tseries = '006';
% tseries = '53';
znum = '002';
% znum = '2';

% cd(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
% path = ['C:\Users\bbscott\Documents\Scott Lab\Raw Data\' sessid];
% path = ['H:\Research\Scott\raw data\' sessid];
% cd(path)

rows = 'V35:X35'; % rows in Session Info to write data
% path = ['C:\Users\bbscott\Documents\Scott Lab\Raw Data\' sessid];
% path=['C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\raw data\' sessid];
% assumes Session Info is in the Dropbox\Bruker directory
cd(path)
AnalyzeSession_tseries(sessid,path,tseries);
% cd(path)
[z_reg]=Zreg_sjk1(sessid,tseries,znum);

% cd('C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');
% cd('C:\Users\sjk\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');
AA = [rmsX, rmsY, z_reg];
xlswrite('Session Info.xlsx',AA,'invivo',rows);
%%
% import processed data
sessid = '20201228-616125';
cd(['C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\processed\mouse\' sessid]);

% now import wheel data
% % 
cd('C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker\raw data')
% % Wheel = readtable('20201213-616125-wheel.xlsx','Range','L2:M812'); 
% Wheel = readtable('wheeldata-20201228macros.xlsm','Sheet','tseries 6','Range','O48:P436'); 
Wheel = readtable('wheeldata-20201228macros.xlsm','Sheet','tseries 4','Range','O69:P470'); 
% % % % 20201213-189
% Wheel = readtable('20201213-616125-wheel.xls','Range','D1:F871'); 
Wheeltime=table2array(Wheel(:,2));
% % 
Wheelspeed= table2array(Wheel(:,1));
%% mouse velocity (x) by total frame displacement (y) scatter plot

xdisp = (x_values-x_values(1)).^2;
ydisp = (y_values-y_values(1)).^2 ; zdisp = (z_values-z_values(1)).^2;
displace = sqrt(xdisp + ydisp + zdisp(1:length(xdisp)));

time = 1:length(displace);duration=67/length(displace);time=time*duration; %sampling rate for resonant galvos is 30fps
velocity = resample(Wheelspeed,5,1);
% resample sometimes puts negative values in areas that should not have
% negative values. manually removed 3 values that were large and then took
% absolute value
velocity=abs(velocity);

vscatter = scatter(velocity(1:length(displace)),displace',30,'MarkerFaceColor','b','MarkerEdgeColor','b');
xlabel('Mouse Velocity in cm/s'); ylabel('Total Focal Plane Displacement (XY) in Micrometers');
alpha(vscatter,.2);
%%

xdisp = (x_values-x_values(1)); ydisp = (y_values-y_values(1)); zdisp = (z_values-z_values(1));
% calculate RMS displacement
    xdispX=xdisp-mean(xdisp);
    ydispY=xdisp-mean(ydisp);
    zdispZ=zdisp-mean(zdisp);
    arithmetic_mean_of_squares=sum(xdispX.^2);
    disprmsX=sqrt(arithmetic_mean_of_squares/length(xdispX));
    arithmetic_mean_of_squares=sum(ydispY.^2);
    disprmsY=sqrt(arithmetic_mean_of_squares/length(ydispY));
    arithmetic_mean_of_squares=sum(zdispZ.^2);
    disprmsZ=sqrt(arithmetic_mean_of_squares/length(zdispZ));
%%
figure(99) 
hold on
xdispscat = scatter(velocity(1:length(xdisp)),xdisp',30,'MarkerFaceColor','r','MarkerEdgeColor','r');
ydispscat = scatter(velocity(1:length(ydisp)),ydisp',30,'MarkerFaceColor','b','MarkerEdgeColor','b');
zdispscat = scatter(velocity(1:length(zdisp)),zdisp',30,'MarkerFaceColor','k','MarkerEdgeColor','k');
% xdispscat = scatter(velocity(1:length(xdisp)),xdisp',30,'MarkerFaceColor','r','MarkerEdgeColor','k');
% ydispscat = scatter(velocity(1:length(ydisp)),ydisp',30,'MarkerFaceColor','b','MarkerEdgeColor','k');
% zdispscat = scatter(velocity(1:length(zdisp)),zdisp',30,'MarkerFaceColor','k','MarkerEdgeColor','k');
xlabel('Mouse Velocity in cm/s'); ylabel('Displacement in Micrometers');
alpha(xdispscat,0.2);alpha(ydispscat,0.2);alpha(zdispscat,0.2);

legenddivert = strcat('Vertical, rms = ' , num2str(disprmsY),'µm');
legenddihoriz = strcat('Horizontal, rms = ',num2str(disprmsX),'µm');
legenddiz=strcat('Z-plane, rms = ',num2str(disprmsZ),'µm');
legend(legenddivert,legenddihoriz,legenddiz);
%% RMS displacement in states of run/rest
% first divide data into periods of run and rest

% threshold
threshold = 0.5; % in cm/s
% xdisp = (x_values-x_values(1)).^2;
% ydisp = (y_values-y_values(1)).^2 ; zdisp = (z_values-z_values(1)).^2;
% displace = sqrt(xdisp + ydisp + zdisp(1:length(xdisp)));
smoothvel = smoothdata(velocity,'gaussian',25);
smoothvel = smoothvel(1:length(displace));
rest = zeros(length(smoothvel),1); run = zeros(length(smoothvel),1);
for ii =1:length(smoothvel) % get two arrays telling you the indices of run and rest
    if (smoothvel(ii) <= 0.15) == 1
        rest(ii) = 1; 
    else
        run(ii) = 1;
    end
end

% now get indices of ones in each rest and run. this tells you the indices 
% to call from smoothvel
krest = find(rest); krun = find(run);

for jj = 1:length(krest) % put velocity values of each index in second column
    krest(jj,2) = smoothvel(krest(jj));
    krest(jj,3) = displace(krest(jj)); % get total displacement values of each index in third column
    krest(jj,4) = xdisp(krest(jj)); % column 4 is x displacement
    krest(jj,5) = ydisp(krest(jj)); % column 5 is y displacement
    krest(jj,6) = zdisp(krest(jj)); % column 6 is z displacement
end
for ll = 1:length(krun)
    krun(ll,2) = smoothvel(krun(ll)); 
    krun(ll,3) = displace(krun(ll));
    krun(ll,4) = xdisp(krun(ll));
    krun(ll,5) = ydisp(krun(ll));
    krun(ll,6) = zdisp(krun(ll));
end

% now calculate RMS displacement
    krunX=krun(:,4)-mean(krun(:,4));
    krunY=krun(:,5)-mean(krun(:,5));
    krunZ=krun(:,6)-mean(krun(:,6));
    arithmetic_mean_of_squares=sum(krunX.^2);
    runrmsX=sqrt(arithmetic_mean_of_squares/length(krunX));
    arithmetic_mean_of_squares=sum(krunY.^2);
    runrmsY=sqrt(arithmetic_mean_of_squares/length(krunY));
    arithmetic_mean_of_squares=sum(krunZ.^2);
    runrmsZ=sqrt(arithmetic_mean_of_squares/length(krunZ));

% now calculate RMS displacement for rest
    krestX=krest(:,4)-mean(krest(:,4));
    krestY=krest(:,5)-mean(krest(:,5));
    krestZ=krest(:,6)-mean(krest(:,6));
    arithmetic_mean_of_squares=sum(krestX.^2);
    restrmsX=sqrt(arithmetic_mean_of_squares/length(krestX));
    arithmetic_mean_of_squares=sum(krestY.^2);
    restrmsY=sqrt(arithmetic_mean_of_squares/length(krestY));
    arithmetic_mean_of_squares=sum(krestZ.^2);
    restrmsZ=sqrt(arithmetic_mean_of_squares/length(krestZ));

% put it in a nice table
RMSdisplaceRestRun = {}; RMSdisplaceRestRun{2,1} = 'Rest';
RMSdisplaceRestRun{3,1} = 'Run';
RMSdisplaceRestRun{1,2} = 'X';RMSdisplaceRestRun{1,3} = 'Y';RMSdisplaceRestRun{1,4} = 'Z';
RMSdisplaceRestRun{2,2} = restrmsX; 
RMSdisplaceRestRun{2,3} = restrmsY; 
RMSdisplaceRestRun{2,4} = restrmsZ; 
RMSdisplaceRestRun{3,2} = runrmsX; 
RMSdisplaceRestRun{3,3} = runrmsY; 
RMSdisplaceRestRun{3,4} = runrmsZ; 


%%
figure12=figure(12);
hold on
sub(1)=subplot(2,1,1);

time = 1:length(x_values);duration=67/length(x_values);time=time*duration; %sampling rate for resonant galvos is 30fps

hold on
smx = smoothdata(x_values,'gaussian',25);
% plot(1:numel(smx),smx,'r','LineWidth',3);
plot(time(1:numel(smx)),smx,'r','LineWidth',3);
smy = smoothdata(y_values,'gaussian',25);
plot(time(1:numel(smx)),smy,'b','LineWidth',3);
smz = smoothdata(z_values,'gaussian',25);
plot(time(1:numel(smx)),smz(1:numel(smx)),'k','LineWidth',3);

xpoints = plot(time(1:numel(x_values)),x_values,'r','LineWidth',3);
xpoints.Color(4)=0.25;
ypoints = plot(time(1:numel(x_values)),y_values,'b','LineWidth',3);
ypoints.Color(4)=0.25;
zpoints=plot(time(1:numel(x_values)),z_values(1:length(x_values)),'k','LineWidth',3);
zpoints.Color(4)=0.25;
% plot(1:numel(x_values),abs(x_values),'LineWidth',3);
% scatterx=scatter(1:numel(x_values),abs(x_values),30,'MarkerFaceColor','b','MarkerEdgeColor','b');
% alpha(scatterx,.2);

ylim([-4 4])
xlim([0 67])
ylabel('Position in micrometers')
% xticklabels({'0','12','24','36','48','60','72','84','96','108','120'})
xlabel('Seconds')
title('in vivo Position Measurements 20201228 t6 z3 Gaussian 25');

legendPosvert = strcat('Vertical, rms = ' , num2str(rmsY),'µm');
legendPoshoriz = strcat('Horizontal, rms = ',num2str(rmsX),'µm');
legendz=strcat('Z-plane, rms = ',num2str(z_reg),'µm');
legend(legendPosvert,legendPoshoriz,legendz);


sub(2)=subplot(2,1,2);

plot(((Wheeltime-Wheeltime(1))/1000),Wheelspeed,'LineWidth',3);
ylabel('Mouse Speed (cm/s)')
xlabel('Seconds')
xlim([0 67])
linkaxes(sub,'x');




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