clear;clc;
addpath(genpath('parser'))
addpath(genpath('animate'))
addpath('quaternions')
%% load inpt files
 [skel,mot] = readMocap('data\HDM05_cut_amc\punchLSide2Reps\HDM_bd.asf','data\HDM05_cut_amc\punchLSide2Reps\HDM_bd_punchLSide2Reps_003_120.amc'); % read data
% standard skeleton coordinate
set(0,'DefaultFigureVisible', 'on')
ratio=2.5/100; % the ratio between the coordinates of data and real world
startframe=1;
endframe=360;  % limit the frame that we want to explore
n1=1;
n2=1;
for f=(startframe:endframe)
    for i=1:24 % reduced to 24 nodes
      index=cell2mat(mot.nameMap(i,3)); % namemap: 
      x(i,n1)=ratio*mot.jointTrajectories{index,1}(1,f);
      y(i,n1)=ratio*mot.jointTrajectories{index,1}(2,f);
      z(i,n1)=ratio*mot.jointTrajectories{index,1}(3,f);
    end
    n1=n1+1;
end
x([15,20],:)=[];  % node 12 15 20 all refer to the same one, so only keep node 12
y([15,20],:)=[];
z([15,20],:)=[];
    % node sequence is shown in an excel
    % parant id mapping
    parentmap=[0,1,2,3,4,1,6,7,8,1,10,11,12,13,12,15,16,17,12,19,20,21];
    
    % ellipsoid parameter 
    a=[0.07,0.07,0.06,0.04,0.07,0.07,0.06,0.04,0.1,0.1,0.1,0.04,0.06,0.05,0.05,0.04,0.04,0.05,0.05,0.04,0.04];
    b=[0.07,0.07,0.06,0.04,0.07,0.07,0.06,0.04,0.1,0.1,0.1,0.04,0.06,0.05,0.05,0.04,0.04,0.05,0.05,0.04,0.04];
for f=(startframe:endframe)
    % build model with ellipsoid
    for k=2:22 
%       plot3([x(k,f),x(parentmap(k),f)], ...
%              [y(k,f),y(parentmap(k),f)], ...
%              [z(k,f),z(parentmap(k),f)])  % connect current node and its parent, so k starts from 2 
    
         current_node(k,:)=[x(k,n2),y(k,n2),z(k,n2)];
         parent_node(k,:)=[x(parentmap(k),n2),y(parentmap(k),n2),z(parentmap(k),n2)];
        c(k-1)=sqrt((x(k,n2)-x(parentmap(k),n2))^2+(y(k,n2)-y(parentmap(k),n2))^2+...
         (z(k,n2)-z(parentmap(k),n2))^2)/2;
         [X,Y,Z]=ellipsoid2P(current_node(k,:),parent_node(k,:),a(k-1), ...  % build the model
             b(k-1),c(k-1),20);
%          surf(X,Y,Z)   
%          
%          hold on
    end
    n2=n2+1;
end

% radar parameter
fs =24e9; % start frequency
Tsweep = 0.375; % Sweep time in ms
Tsweep=Tsweep/1000; %then in sec
NTS =256 ; % Number of time samples per sweep
nsequence=[1:NTS];
Fs=NTS./Tsweep;
Bw = 4*10^9; % FMCW Bandwidth
S=Bw./Tsweep;  % slope
radarloc=[-4,1,1]; % radar location
c_speed=3e8; % light speed
lamda=c_speed/fs; % wavelength
% interpolation among nodes
inx=[1:endframe-startframe+1];
numchirp=floor((endframe-startframe+1)./(120*Tsweep)); % framerate=120hz, number of chirps
totalsample=NTS*numchirp;
xp=[1:(endframe-startframe)/(totalsample-1):endframe-startframe+1];  % interpolation position
for k=(1:22) 
    intrp_x(k,:)=pchip(inx,x(k,:),xp);
    node(k).x=reshape(intrp_x(k,:),NTS,numchirp); % interpolate and reshape (NTS,numchirp)
    intrp_y(k,:)=pchip(inx,y(k,:),xp);
    node(k).y=reshape(intrp_y(k,:),NTS,numchirp);
    intrp_z(k,:)=pchip(inx,z(k,:),xp);
    node(k).z=reshape(intrp_z(k,:),NTS,numchirp);
end
addwindow1=hanning(NTS);
addwindow2=hanning(numchirp);
% test the shape after inteprolation

% for f=3:3
%     for k=2:22
%         plot3([intrp_x(k,f),intrp_x(parentmap(k),f)], ...
%          [intrp_y(k,f),intrp_y(parentmap(k),f)], ...
%          [intrp_z(k,f),intrp_z(parentmap(k),f)])  % connect current node and its parent
%         
%         current_node(k,:)=[intrp_x(k,f),intrp_y(k,f),intrp_z(k,f)];
%         parent_node(k,:)=[intrp_x(parentmap(k),f),intrp_y(parentmap(k),f),intrp_z(parentmap(k),f)];
%         c(k-1)=sqrt((intrp_x(k,f)-intrp_x(parentmap(k),f))^2+(intrp_y(k,f)-intrp_y(parentmap(k),f))^2+...
%         (intrp_z(k,f)-intrp_z(parentmap(k),f))^2)/2;
%         [X,Y,Z]=ellipsoid2P(current_node(k,:),parent_node(k,:),a(k-1), ...
%          b(k-1),c(k-1),20);
%         surf(X,Y,Z)
%         hold on
%     end
% end
for k=(2:22) % node index
    for nc=(1:numchirp) % chirp index
            for ns=(1:NTS)  % sample index
            middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
            middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
            middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
            distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
            % aspect vector of the current node
            aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                   node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
            % calculate theta angle
            A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
            B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
            A_dot_B=dot(A,B,1);
            A_sum_sqrt=sqrt(sum(A.*A,1));
            B_sum_sqrt=sqrt(sum(B.*B,1));
            theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
            phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
            rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
            amp= sqrt(rcs); 
            % IF singnal
            phi_IF=4*pi*distance(ns)*fs/c_speed;
            omega_IF=4*pi*S*distance(ns)/(c_speed*Fs);
            skel_radareceive(k-1).data(ns,nc)=amp*exp(j*(phi_IF+omega_IF*(ns-1)));  % received IF signal 
            end
    end
    skel_radareceive(k-1).data= awgn(skel_radareceive(k-1).data,20);
end


% test distance change between 1st and last frame of each skeleton
% figure
% for nc=[1,numchirp]
%     inty=abs(fft(skel_radareceive(12).data(:,nc),256));
%     xaxis=[1:128];
%     r=xaxis*c_speed*Fs/(2*S*NTS);
%     plot(r,inty)
%     xlabel("Distance(m)")
%     ylabel("Intensity")
%     title("neck-head")
%     legend("The first chirp", "The last chirp")
%     hold on
% end

% range profile
time_axis=[1:numchirp];
range_axis=nsequence*c_speed*Fs/(2*S*NTS);
data_power=zeros(NTS,numchirp);
data=zeros(NTS,numchirp);
for nc=(1:numchirp)
    for k=(1:21)
        skel_radareceive(k).data(:,nc)= skel_radareceive(k).data(:,nc).*addwindow1;
        data_power(:,nc)=abs(fft(skel_radareceive(k).data(:,nc),256))+data_power(:,nc);  % for range calculation 
        data(:,nc)=fft(skel_radareceive(k).data(:,nc),256)+data(:,nc);  % for velocity calculation
    end
end

figure
colormap(jet(256))
imagesc(time_axis*Tsweep,range_axis,data_power)
axis xy
xlabel('Time(s)','FontSize',16)
ylabel('Range(m)','FontSize',16) 
%title("Combination")
colorbar 
drawnow
hold on

% vmax=lamda/4/Tsweep;
% doppler_axis = linspace(vmax,-vmax,numchirp);
% % This selects the range bins where we want to calculate the spectrogram
% bin_indl = 10;
% bin_indu = 30; % 128 samples of each chirp, only take 10-100bins
% 
% micro_doppler.PRF=1/Tsweep;
% micro_doppler.TimeWindowLength=128;
% micro_doppler.OverlapFactor = 0.95;
% micro_doppler.OverlapLength = round(micro_doppler.TimeWindowLength*micro_doppler.OverlapFactor);
% micro_doppler.Pad_Factor = 1;
% micro_doppler.FFTPoints = micro_doppler.Pad_Factor*micro_doppler.TimeWindowLength;
% micro_doppler.DopplerBin=micro_doppler.PRF/(micro_doppler.FFTPoints);
% micro_doppler.DopplerAxis=-micro_doppler.PRF/2:micro_doppler.DopplerBin:micro_doppler.PRF/2-micro_doppler.DopplerBin;
% micro_doppler.WholeDuration=size(data,2)/micro_doppler.PRF;
% Data_spec_MTI2=0;
% 
% % STFT to extract micro-doppler information
% 
% % for RBin=bin_indl:1:bin_indu
% %     Data_MTI_temp = fftshift(spectrogram(data(RBin,:),micro_doppler.TimeWindowLength,micro_doppler.OverlapLength,micro_doppler.FFTPoints),1);
% %     Data_spec_MTI2=Data_spec_MTI2+abs(Data_MTI_temp);                                
% % end
% % micro_doppler.TimeAxis=linspace(0,micro_doppler.WholeDuration,size(Data_spec_MTI2,2));
% % 
% % Data_spec_MTI2=flipud(Data_spec_MTI2); % due to fftshift, doppler axis is from negative to positive
% % % Data_spec_MTI2((length(doppler_axis)+1)/2,:)=0;
% % figure %show the doppler spectrum
% % imagesc(micro_doppler.TimeAxis,micro_doppler.DopplerAxis.*3e8/2/fs,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
% % % ylim([-6 6]); colorbar
% % colormap; %xlim([1 9])
% % clim = get(gca,'CLim');
% % set(gca, 'CLim', clim(2)+[-40,0]);
% % xlabel('Time[s]', 'FontSize',16);
% % ylabel('Velocity [m/s]','FontSize',16)
% % %title("HDM_dg_06-03_03-120","Interpreter","none")
% % set(gca, 'FontSize',16)
% 
% % Range-doppler matrix doppler DFT
% for n=(1:NTS)
%     data(n,:)=data(n,:).*addwindow2';
%     RDMabs(n,:)=abs(fftshift(fft(data(n,:),numchirp)));
%     RDMcom(n,:)=(fftshift(fft(data(n,:),numchirp)));
% end
% % 
% RDM=20*log10(RDMabs);
% r = 7;    % 滤波半径
% a = 2;    % 全局方差
% b = 1;  % 局部方差
% 
% [~,~,ch] = size(RDM);
% % 判断是灰度图还是彩色图像
% 
% % RDM = bfilt_gray(RDM,r,a,b);
% 
% figure
% imagesc(doppler_axis,range_axis,RDM); colormap('jet'); 
% % xlim([-6 6]); colorbar
% xlim([-vmax vmax]); colorbar
% ylim=[0 10];
% colormap; %xlim([1 9])
% clim = get(gca,'CLim');
% set(gca, 'CLim', clim(2)+[-40,0]);
% xlabel('Velocity(m/s)','FontSize',16);
% ylabel('Range(m)','FontSize',16)
% axis xy %the y-axis is vertical with values increasing from bottom to top.
% set(gca, 'FontSize',16)
% 
% % % 2d_ca_cfar and display the CFAR output using the Surf function 
% 
% [CUT,det_rangeindex,det_veloindex]=ca_cfar(RDM); %这里的rangeindex是从上往下数的，veloindex是从左往右
% figure('Name', 'CA-CFAR Filtered RDM')
% surf(doppler_axis,range_axis,CUT);
% xlim([-vmax vmax]); 
% ylim=[0 10];
% colorbar;
% %title( 'CA-CFAR Filtered RDM surface plot');
% xlabel('Speed(m/s)');
% ylabel('Range(m)');
% zlabel('Normalized Amplitude');
% set(gca, 'FontSize',16)
% 
% figure('Name', 'CA-CFAR Filtered Range')
% surf(doppler_axis,range_axis,CUT);
% colorbar;
% title( 'CA-CFAR Filtered RDM surface plot');
% xlabel('Speed');
% ylabel('Range');
% zlabel('Normalized Amplitude');
% view(90,0);
% 
% figure('Name', 'CA-CFAR Filtered Speed')
% surf(doppler_axis,range_axis,CUT);
% colorbar;
% title( 'CA-CFAR Filtered RDM surface plot');
% xlabel('Speed');
% ylabel('Range');
% zlabel('Normalized Amplitude');
% view(0,0);
% % 
% % 2d_os_cfar and display the CFAR output using the Surf function 
% % [CUT,det_rangeindex,det_veloindex]=os_cfar(RDM);
% % figure('Name', 'OS-CFAR Filtered RDM')
% % surf(doppler_axis,range_axis,CUT);
% % xlim([-vmax vmax]); 
% % ylim=[0 10];
% % colorbar;
% % title( 'OS-CFAR Filtered RDM surface plot');
% % xlabel('Speed(m/s)');
% % ylabel('Range(m)');
% % zlabel('Normalized Amplitude');
% % set(gca, 'FontSize',16)
% % 
% % figure('Name', 'OS-CFAR Filtered Range')
% % surf(doppler_axis,range_axis,CUT);
% % colorbar;
% % title( 'OS-CFAR Filtered RDM surface plot');
% % xlabel('Speed(m/s)');
% % ylabel('Range(m)');
% % zlabel('Normalized Amplitude');
% % view(90,0);
% % 
% % figure('Name', 'OS-CFAR Filtered Speed')
% % surf(doppler_axis,range_axis,CUT);
% % colorbar;
% % title( 'OS-CFAR Filtered RDM surface plot');
% % xlabel('Speed(m/s)');
% % ylabel('Range(m)');
% % zlabel('Normalized Amplitude');
% % view(0,0);
% target=0;
% noise=0;
% dis=zeros(size(RDM,1),size(RDM,2));
% for j=1:size(RDM,1)
%     for k=1:size(RDM,2)
%         if RDM(j,k)>=50
%             RDM(j,k)=1;
%             target=target+1;
%         else
%             RDM(j,k)=0;
%             noise=noise+1;
%         end
%         dis(j,k)=CUT(j,k)-RDM(j,k);
%     end
% end
% 
% fn=length(find(dis==-1));
% fp=length(find(dis==1));
% tp=target-fn;
% presision=tp./(tp+fp);
% fnr=fn./target; 
% fpr=fp./noise;
% disp(presision)
% disp(fnr)