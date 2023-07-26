clear;clc;
addpath(genpath('parser'))
addpath(genpath('animate'))
addpath('quaternions')
delete(gcp('nocreate'));
% load inpt files
[skel,mot] = readMocap('data\HDM05_cut_amc\jumpingJack3Reps\HDM_bd.asf','data\HDM05_cut_amc\jumpingJack3Reps\HDM_bd_jumpingJack3Reps_001_120.amc'); % read data
set(0,'DefaultFigureVisible', 'on');
actyclass='jumpingJack3Reps';
tic   
parpool(8);
framerate=120;
window_size=12; %设置窗口大小及滑动长度
sliding_size=3;
period=120;
n=0;
frame_limit= floor((period-window_size)/sliding_size);
time=((0:frame_limit)*sliding_size+period*n)/framerate;
test_start_frame=1+period*n;
test_end_frame=window_size+period*n;
[range,velocity,CUT]=RDS_method(test_start_frame,test_end_frame,mot);
% [X,Y]=meshgrid(velocity,range);
% contour(X,Y,70*CUT,[70,70])
% m=contour(X,Y,70*CUT,[70,70]);  % multiple 70 to avoid range equal to the intensity
% for k =(1:size(m,2))
%     if(m(1,k)==70)
%        continue
%     else
%         k1=find(range==m(1,k));
%         k2=find(velocity==m(2,k));
%         new_rdm(k1,k2)=1;
%     end
% end                     提取轮廓
V(:,1,:)=CUT;

parfor i =1:frame_limit
   startframe=period*n+1+i*sliding_size;
   endframe=period*n+window_size+i*sliding_size;
   [~,~,CUT]=RDS_method(startframe,endframe,mot);
%    [X,Y]=meshgrid(velocity,range);
%     contour(X,Y,70*CUT,[1,1])
%     t=contour(X,Y,70*CUT,[1,1]);
%     for k =(1:size(t,2))
%         if(t(1,k)==70)
%             continue
%         else
%             m1=find(range==t(1,k));
%             m2=find(velocity==t(2,k));
%             new_rdm(m1,m2)=1;
%         end
%     end
    V(:,i+1,:)=CUT;
end
%figure
iso_value=0.5;
[nx,ny,nz]=meshgrid(time,range,velocity); % return size y-x-z (range,time,velocity)
%isosurface(nx,ny,nz,V,iso_value)
camlight;
lighting gouraud;
[f,v]=isosurface(nx,ny,nz,V,iso_value);
num_of_face=size(f,1);
actyclass='punchLSide2Reps';
for t=1:2
    K=randsample(num_of_face,3000);
    center=zeros(length(K),3);
    filename=strcat('actyclass\',actyclass,'\',datestr(now,30),'.pts');
    fid = fopen(filename, 'w');
    for i=1:length(K)
        center(i,:)=(v(f(K(i),1),:)+v(f(K(i),2),:)+v(f(K(i),3),:))/3;
        fprintf(fid,'%f %f %f\n', center(i,1),center(i,2),center(i,3));
    end
    fclose(fid);
end

% figure
% scatter3(center(:,1),center(:,2),center(:,3),'filled')
% xlabel('Time(s)')
% ylabel('Range(m)')
% zlabel('Velocity(m/s)')
% toc



