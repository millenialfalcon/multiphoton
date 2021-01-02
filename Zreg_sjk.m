function [z_reg]=Zreg_sjk1(sessid,tseries,znum)
%clear; getsessids; sessid=sessids{3};
close all
%kinematic
% path='C:\Users\bbscott\Documents\Scott Lab\processed data\';
% cd([path sessid])
% load(['RegData' sessid '-' tseries '.mat'])
ppath = ['C:\Users\bbscott\Documents\Scott Lab\processed data\' sessid];
% path = ['C:\Users\bbscott\Documents\Scott Lab\Raw Data\' sessid];
cd(ppath)



% path = ['H:\Research\Scott\processed data\headplate\mouse\' sessid];
% path = 'C:\Users\bbscott\Downloads\' sessid];
% cd(path)

load(['RegData' sessid '-' tseries '.mat'])
% path = ['H:\Research\Scott\raw data\' sessid];
path = ['C:\Users\bbscott\Documents\Scott Lab\Raw Data\' sessid];
cd(path)


% path = ['H:\Research\Scott\raw data\' sessid];
% cd(path)

Y=[]; Z=[];Zpos=[];%path=pwd;
%% first open the Zseries
folder=dir(['ZSeries*' znum]);
cd(folder(1).name)
% files=dir('*0511*');
files=dir('*Ch1*'); % list files only with ch2 in name % these are the files you want
ch1files=dir('*Ch2*');
for k=1:length(files)
    Z1 = read_file(files(k).name);
    if k==1
        data = bfopen(files(k).name); %Read the header
        omeMeta = data{1, 4};%Extract metadata
        Z=Z1;
    else
        Z=cat(3,Z,Z1);
    end
    if isempty(ch1files)
        planePosZ =  omeMeta.getPlanePositionZ(0,k-1); %compute Z position
    elseif ~isempty(ch1files)
        planePosZ =  omeMeta.getPlanePositionZ(0,2*k-1);
    end
    Zpos=[Zpos planePosZ.value()];
    
end
Z = single(Z); % convert to single precision
T = size(Z,ndims(Z));
Z = Z - min(Z(:));

%% Compute Z locations
mpp=omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
microns_per_pixel=double(mpp);
maxshift=round(10/microns_per_pixel);
options_rigid = NoRMCorreSetParms('d1',size(Z,1),'d2',size(Z,2),'bin_width',200,'max_shift',maxshift,'us_fac',50,'init_batch',200);
mov=M1; %here is where you select the movie!  M1 is the registered movie
thisshiftgauss=nan(1,size(mov,3));
Rs=nan(size(Z1,3)-10,size(mov,3));
gaussEqn = 'a*exp(-((x-b)/c)^2)+d'; 
% startpoints=[1 0 1 0.1];
startpoints=[0.6 -0.3 0.6 -0.3];
% startpoints=startpoints; % if f1 returns 0, try different
                         % startpoints
template=squeeze(mean(mov,3));
[Z1,~,~,~] = normcorre(Z,options_rigid,template);
for j=1:size(mov,3)
    for i=1:(size(Z1,3))
        R = corrcoef(mov(:,:,j),Z1(:,:,i)); % calculate correlation coefficient 
        Rs(i,j)=R(1,2);
    end
    
    
    y=Rs(:,j);
    x=double(Zpos);
    f1 = fit(x,y,gaussEqn,'Start', startpoints');
    
 if mean(Rs(:,j))>0.1  %we exclude points with low correlation values.
    thisshiftgauss(j)=f1.b;
%     figure(1); plot(f1,x,y)
%         title(num2str(f1.b));
%     ylim([0.3 0.6]);
%     drawnow
 else
     thisshiftgauss(j)=nan;
 end
end
%%
z_values=thisshiftgauss-thisshiftgauss(1);
figure;
x=double(Zpos);
figure;
subplot(1,2,1)
plot(z_values,'.','Markersize',20)
ylabel('Z position')
xlabel('Minutes')

subplot(1,2,2)
imagesc(flipud(Rs))
ylabel('Zposition')
xlabel('Minutes')

figure(3); scatter(x_values,y_values)
ylabel('Vertical displacement (um)');xlabel('Horizontal displacement (um)');


X=thisshiftgauss-nanmean(thisshiftgauss); %error
arithmetic_mean_of_squares=nansum(X.^2);
errorinZ=sqrt(arithmetic_mean_of_squares/length(X));
z_reg=sqrt(arithmetic_mean_of_squares/length(X));
% cd(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
cd(['C:\Users\bbscott\Documents\Scott Lab\processed data\' sessid])
save(['RegData' sessid '-' tseries '-' znum '.mat'],'M1','x_values','y_values',...
'rmsX','rmsY','microns_per_pixel','z_values','thisshiftgauss','z_reg','errorinZ','-v7.3')


