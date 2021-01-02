function AnalyzeSession_tseries(sessid, path,tseries)
%Compute registration accuracy
cd(path)
%   cd(['TSeries-12192020-1101-001' sessid])

Y=[];Xpos=[];Ypos=[]; 
y_values=[]; x_values=[];
folders=dir(['TSeries*' tseries]);
cd(folders(1).name); %go to the first folder
% files=dir('*189*');
files=dir('*Ch1*'); % list all of the files in this folder that have Ch2 in the name
    for k=1:length(files)
        Y1 = read_file(files(k).name);
        if k==1
            data = bfopen(files(k).name); %Read the header
        omeMeta = data{1, 4};%Extract metadata
            Y=Y1;
        else
            Y=cat(3,Y,Y1);
        end
    end
    Y = single(Y); % convert to single precision
T = size(Y,ndims(Y));
Y = Y - min(Y(:));
mpp=omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER); % get units
microns_per_pixel=double(mpp); %datatype shit
maxshift=round(15/microns_per_pixel);
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',maxshift,'us_fac',50,'init_batch',200);

% perform motion correction
template=Y(:,:,1); %template is the average of the first 10 frames
[M1,shifts1,template1,~] = normcorre(Y,options_rigid,template);
M0=Y;
for i=1:size(M1,3)
    R = corrcoef(template1,M1(:,:,i));
    Rs(i)=R(1,2);
end
mpp=omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROMETER);
imagetime=omeMeta.getImageAcquisitionDate(0); % tells you the time on X axis
numvals=length(shifts1);

for i=1:numvals
    y_values(i)=shifts1(i).shifts(:,:,1,1).*microns_per_pixel;
    x_values(i)=shifts1(i).shifts(:,:,1,2).*microns_per_pixel;
end
clear M0; 
%remove unregistered 
y_values=y_values(Rs>(mean(Rs(2:end))-2*std(Rs(2:end))));
x_values=x_values(Rs>(mean(Rs(2:end))-2*std(Rs(2:end))));

subplot(1,2,2); hold on
plot(x_values,'.r')
plot(y_values,'.b')
ylim([-7 7])
xlabel('Frames')
ylabel('Position (microns)')
title(['Reg for ' sessid]);
% compute RMS displacement
X=x_values-mean(x_values);
Y=y_values-mean(y_values);
arithmetic_mean_of_squares=sum(X.^2);
rmsX=sqrt(arithmetic_mean_of_squares/length(x_values));
arithmetic_mean_of_squares=sum(Y.^2);
rmsY=sqrt(arithmetic_mean_of_squares/length(y_values));

legend(['rmsX=' num2str(rmsX)],['rmsY=' num2str(rmsY)], 'Location','northeast')
%% save the data
fr=mean(Rs<(mean(Rs(2:end))-2*std(Rs(2:end))));
fr=round(fr,2);


if ~exist(['C:\Users\bbscott\Documents\Scott Lab\processed data\' sessid],'dir')
    mkdir(['C:\Users\bbscott\Documents\Scott Lab\processed data\' sessid])
end
    cd(['C:\Users\bbscott\Documents\Scott Lab\processed data\' sessid])
    
    
% if ~exist(['H:\Research\Scott\processed data\headplate\mouse\' sessid],'dir')
%     mkdir(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
% end
%     cd(['H:\Research\Scott\processed data\headplate\mouse\' sessid])
%     
    
% saveas(gcf,['Data for' sessid '.pdf'])
% save(['RegData' sessid '.mat'],'M0','M1','D0','x_values','y_values')
save(['RegData' sessid '-' tseries '.mat'],'imagetime','M1','x_values','y_values',...
    'rmsX','rmsY','microns_per_pixel','-v7.3') 
% errorinz=computeZ(sessid); 
% z_drift=Zdrift_sjk(sessid); 
% z_reg=Zreg_sjk1(sessid,tseries,znum);
% cd('C:\Users\bbscott\Dropbox (BOSTON UNIVERSITY)\Kinematic_Data\Bruker');
% AA = [rmsX, rmsY, z_reg];
% xlswrite('Session Info.xlsx',AA,'Sheet1',rows);
%need to use v7.3 for variables >2GB (M0,M1)
