clear,clc;
delete(gcp('nocreate'));
set(0,'DefaultFigureVisible', 'off')
n=parpool(6);
addpath(genpath('parser'))
addpath(genpath('animate'))
addpath('quaternions')
subject="dg";
actyclass="walkLeft3Steps";
sliding_window=12;
tic
[skel,mot] = readMocap('data\HDM05_cut_amc\walkLeft3Steps\HDM_dg.asf','data\HDM05_cut_amc\walkLeft3Steps\HDM_dg_walkLeft3Steps_008_120.amc'); % read data
%[skel,mot] = readMocap('data/HDM_dg.asf', 'data/HDM_dg_06-03_03_120.amc');
max=floor(mot.nframes/sliding_window)-1;
parfor i=(0:max)
    startframe=1+i*sliding_window;
    endframe=12+i*sliding_window; % 0.1 seconds for each radar frame
    cloudpoints=fdradar(startframe,endframe,mot); 
    filename=strcat('acty_add_velocity\',actyclass,'\',subject,'_',int2str(8),'_',int2str(i),'.pts');
    mat_to_pts(cloudpoints,filename)
%     figure
%     disp(size(cloudpoints.x))
%     scatter3(cloudpoints.x,cloudpoints.y,cloudpoints.z,10,'filled')
%     hold on 
end
toc