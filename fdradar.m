function cloudpoints=fdradar(startframe,endframe,mot)
    % standard skeleton coordinate
    set(0,'DefaultFigureVisible', 'off')
    ratio=2.5/100; % the ratio between the coordinates of data and real world
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
    fs =12e9; % start frequency
    Tsweep = 0.0625 ; % Sweep time in ms
    Tsweep=Tsweep/1000; %then in sec
    NTS =256 ; % Number of time samples per sweep
    nsequence=[1:NTS];
    Fs=NTS./Tsweep;
    Bw = 4e9; % FMCW Bandwidth
    S=Bw./Tsweep;  % slope
    radarloc=[-4,1,1]; % radar location
    c_speed=3e8; % light speed
    lamda=c_speed/fs; % wavelength
    numty=6;
    numtz=6;% number of transmitters
    numry=8;
    numrz=8;% number of receivers
    numt=numty+numtz;
    dy=lamda/2;
    dz=lamda/2;
    vitualarray=numty*numry+numtz*numrz;
    
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
            if mod(nc,numt)==1 
                ty=0; % tx0
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+11)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==2
                ty=1; %tx1
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))./sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+10)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==3
                    ty=2; % tx0
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+9)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
             elseif mod(nc,numt)==4
                    ty=3; % tx0
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+8)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==5
                    ty=4; % tx0
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(-2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+7)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==6
                    ty=5; % tx0
                for ry=(0:numry-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_azi=(middlepoint(k).y(ns,nc)-radarloc(2))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(ty*numry+ry)*dy*sin_azi);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+6)/numt,ty*numry+ry+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==7
                tz=0;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                   sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+5)/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end   
            elseif mod(nc,numt)==8
                tz=1;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*j*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*j*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*j/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+4)/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==9
                tz=2;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*j*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*j*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*j/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+3)/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==10
                tz=3;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+2)/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==11
                tz=4;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,(nc+1)/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end
            elseif mod(nc,numt)==0
                tz=5;
                for rz=(0:numrz-1)
                    for ns=(1:NTS)  % sample index
                    middlepoint(k).x(ns,nc)=0.5*[node(parentmap(k)).x(ns,nc)+node(k).x(ns,nc)];
                    middlepoint(k).y(ns,nc)=0.5*[node(parentmap(k)).y(ns,nc)+node(k).y(ns,nc)];
                    middlepoint(k).z(ns,nc)=0.5*[node(parentmap(k)).z(ns,nc)+node(k).z(ns,nc)];
                    distance(ns)=sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+(middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % aspect vector of the current node
                    aspect(:,k)=[node(parentmap(k)).x(ns,nc)-node(k).x(ns,nc); node(parentmap(k)).y(ns,nc)-node(k).y(ns,nc);...
                           node(parentmap(k)).z(ns,nc)-node(k).z(ns,nc)];
                    % calculate rcs
                    A=[radarloc(1)-node(k).x(ns,nc);radarloc(2)-node(k).y(ns,nc);radarloc(3)-node(k).z(ns,nc)];
                    B=[aspect(1,k);aspect(2,k);aspect(3,k)];    
                    A_dot_B=dot(A,B,1);
                    A_sum_sqrt=sqrt(sum(A.*A,1));
                    B_sum_sqrt=sqrt(sum(B.*B,1));
                    theta=acos(A_dot_B./(A_sum_sqrt.*B_sum_sqrt));
                    phi=asin((radarloc(2)-node(k).y(ns,nc))./sqrt((node(k).x(ns,nc)-radarloc(1))^2+(node(k).y(ns,nc)-radarloc(2))^2));
                    rcs=rcs_ellipsoid(a(k-1),b(k-1),c(k-1),phi,theta); % calculate rcs
                    amp= sqrt(rcs); 
                    % calculate the sin(angle)
                    sin_ele=(middlepoint(k).z(ns,nc)-radarloc(3))/sqrt((middlepoint(k).x(ns,nc)-radarloc(1))^2+(middlepoint(k).y(ns,nc)-radarloc(2))^2+ ...
                        (middlepoint(k).z(ns,nc)-radarloc(3))^2);
                    % transmit the signal
                    t0=(ns-1)*1/Fs;
                    st=exp(2*pi*1i*(fs*t0+S*t0^2/2));
                    % receive the signal
                    tau=2*distance(ns)/c_speed;
                    sr=amp*exp(2*pi*1i*(fs*(t0-tau)+S*(t0-tau)^2/2))*exp(2*pi*1i/lamda*(tz*numrz+rz)*dz*sin_ele);
                    % calculate the IF signal
                    skel_radareceive(k-1).data(ns,nc/numt,tz*numrz+rz+numry*numty+1)=st*conj(sr);
                    end
                end
            end
        end
       skel_radareceive(k-1).data= awgn(skel_radareceive(k-1).data,20);
    end
    
    % range profile
    numchirp=floor(numchirp/numt);
    time_axis=2*Tsweep*[1:numchirp];
    range_axis=nsequence*c_speed*Fs/(2*S*NTS);
    data_power=zeros(NTS,numchirp,vitualarray);
    data=zeros(NTS,numchirp,vitualarray);
    addwindow1=hanning(NTS);
    addwindow2=hanning(numchirp);
    for va=(1:vitualarray)
        for k=(1:21)
            for nc=(1:numchirp)
                skel_radareceive(k).data(:,nc,va)= skel_radareceive(k).data(:,nc,va).*addwindow1;
                data_power(:,nc,va)=abs(fft(skel_radareceive(k).data(:,nc,va),NTS))+data_power(:,nc,va);  % for range calculation 
                data(:,nc,va)=fft(skel_radareceive(k).data(:,nc,va),NTS)+data(:,nc,va);  % for velocity calculation
            end
        end
    end
    
    % test distance change between 1st and last frame of each skeleton
    % figure
    % for nc=(1:10)
    %     inty=abs(fft(skel_radareceive(12).data(:,nc,1),128));
    %     xaxis=[1:128];
    %     r=xaxis*c_speed*Fs/(2*S*NTS);
    %     plot(r,inty)
    %     xlabel("Distance(m)")
    %     ylabel("Intensity")
    %     title("neck-head")
    %     legend("The first chirp", "The last chirp")
    %     hold on
    % end
    
    % range detection 
    figure
    colormap(jet(256))
    imagesc(time_axis,range_axis,data_power(:,:,1)) % range detection of vitualarray 1
    axis xy
    xlabel('Number of Chirp')
    ylabel('Distance(m)')
    title("Combination")
    colorbar 
    drawnow
    hold on
    
    vmax=lamda/4/(numt*Tsweep);
    doppler_axis = linspace(vmax,-vmax,numchirp);
    % This selects the range bins where we want to calculate the spectrogram
    bin_indl = 1;
    bin_indu = 128; % 128 samples of each chirp, only take 10-100bins
    
    micro_doppler.PRF=1/Tsweep;
    micro_doppler.TimeWindowLength=128;
    micro_doppler.OverlapFactor = 0.95;
    micro_doppler.OverlapLength = round(micro_doppler.TimeWindowLength*micro_doppler.OverlapFactor);
    micro_doppler.Pad_Factor = 1;
    micro_doppler.FFTPoints = micro_doppler.Pad_Factor*micro_doppler.TimeWindowLength;
    micro_doppler.DopplerBin=micro_doppler.PRF/(micro_doppler.FFTPoints);
    micro_doppler.DopplerAxis=-micro_doppler.PRF/2:micro_doppler.DopplerBin:micro_doppler.PRF/2-micro_doppler.DopplerBin;
    micro_doppler.WholeDuration=size(data,2)/micro_doppler.PRF;
    Data_spec_MTI2=0;
    
    % STFT to extract micro-doppler information
    
    % for RBin=bin_indl:1:bin_indu
    %     Data_MTI_temp = fftshift(spectrogram(data(RBin,:),micro_doppler.TimeWindowLength,micro_doppler.OverlapLength,micro_doppler.FFTPoints),1);
    %     Data_spec_MTI2=Data_spec_MTI2+abs(Data_MTI_temp);                                
    % end
    % micro_doppler.TimeAxis=linspace(0,micro_doppler.WholeDuration,size(Data_spec_MTI2,2));
    % 
    % Data_spec_MTI2=flipud(Data_spec_MTI2); % due to fftshift, doppler axis is from negative to positive
    % figure %show the doppler spectrum
    % imagesc(micro_doppler.TimeAxis,micro_doppler.DopplerAxis.*3e8/2/fc,20*log10(abs(Data_spec_MTI2))); colormap('jet'); axis xy
    % ylim([-6 6]); colorbar
    % colormap; %xlim([1 9])
    % clim = get(gca,'CLim');
    % set(gca, 'CLim', clim(2)+[-40,0]);
    % xlabel('Time[s]', 'FontSize',16);
    % ylabel('Velocity [m/s]','FontSize',16)
    % title("HDM_dg_06-03_03-120","Interpreter","none")
    % set(gca, 'FontSize',16)
    
    % Range-doppler matrix doppler DFT
    for va=(1:vitualarray)
        for n=(1:NTS)
        data(n,:,va)=data(n,:,va).*addwindow2';
        RDMabs(n,:,va)=abs(fftshift(fft(data(n,:,va),numchirp)));
        RDMcom(n,:,va)=(fftshift(fft(data(n,:,va),numchirp)));
        end
    end
    RDM=sum(RDMabs,3);
    RDM=20*log10(RDM);
    r = 7;    % 滤波半径
    a = 2;    % 全局方差
    b = 1;  % 局部方差
    
    [~,~,ch] = size(RDM);
    % 判断是灰度图还是彩色图像
    
    %RDM = bfilt_gray(RDM,r,a,b);
    
    figure
    imagesc(doppler_axis,range_axis,RDM); colormap('jet'); axis xy
    xlim([-vmax vmax]); colorbar
    ylim=[0 10];
    colormap; %xlim([1 9])
    clim = get(gca,'CLim');
    set(gca, 'CLim', clim(2)+[-40,0]);
    xlabel('Velocity(m/s)','FontSize',16);
    ylabel('Range(m)','FontSize',16)
    set(gca, 'FontSize',16)
%     
%     2d_ca_cfar and display the CFAR output using the Surf function 
%     
%     [CUT,det_rangeindex,det_veloindex]=ca_cfar(RDM);
%     figure('Name', 'CA-CFAR Filtered RDM')
%     surf(doppler_axis,range_axis,CUT);
%     colorbar;
%     title( 'CA-CFAR Filtered RDM surface plot');
%     xlabel('Speed(m/s)');
%     ylabel('Range(m)');
%     zlabel('Normalized Amplitude');
%     
%     figure('Name', 'CA-CFAR Filtered Range')
%     surf(doppler_axis,range_axis,CUT);
%     colorbar;
%     title( 'CA-CFAR Filtered RDM surface plot');
%     xlabel('Speed(m/s)');
%     ylabel('Range(m/s)');
%     zlabel('Normalized Amplitude');
%     view(90,0);
%     
%     figure('Name', 'CA-CFAR Filtered Speed')
%     surf(doppler_axis,range_axis,CUT);
%     colorbar;
%     title( 'CA-CFAR Filtered RDM surface plot');
%     xlabel('Speed(m/s)');
%     ylabel('Range(m)');
%     zlabel('Normalized Amplitude');
%     view(0,0);
    
    % 2d_os_cfar and display the CFAR output using the Surf function 
    [CUT,det_rangeindex,det_veloindex]=os_cfar(RDM);
    figure('Name', 'OS-CFAR Filtered RDM')
    surf(doppler_axis,range_axis,CUT);
    colorbar;
    title( 'OS-CFAR Filtered RDM surface plot');
    xlabel('Speed(m/s)');
    ylabel('Range(m)');
    zlabel('Normalized Amplitude');
    
    figure('Name', 'OS-CFAR Filtered Range')
    surf(doppler_axis,range_axis,CUT);
    colorbar;
    title( 'OS-CFAR Filtered RDM surface plot');
    xlabel('Speed(m/s)');
    ylabel('Range(m)');
    zlabel('Normalized Amplitude');
    view(90,0);
    
    figure('Name', 'OS-CFAR Filtered Speed')
    surf(doppler_axis,range_axis,CUT);
    colorbar;
    title( 'OS-CFAR Filtered RDM surface plot');
    xlabel('Speed(m/s)');
    ylabel('Range(m)');
    zlabel('Normalized Amplitude');
    view(0,0);
%     
    
    % phase compensation
    vitualarray_yz=48;
    n=[-vitualarray_yz/2:1:vitualarray_yz/2-1];
    sin_angle_axis=n*lamda./(vitualarray_yz*dy);
    index=[det_rangeindex;det_veloindex];
    % for vx=(5:8)
    %     for col=(-numchirp/2:1:numchirp/2-1)
    %         RDMcom(:,col+numchirp/2+1,vx)=RDMcom(:,col+numchirp/2+1,vx)*exp(-1j*2*pi*col/(numchirp*numtx));
    %     end
    % end
    for col=1:numchirp
        for k=1:11
            RDMcom(:,col,[8*k+1:8*k+8])=RDMcom(:,col,[8*k+1:8*k+8])*exp(-1j*2*pi*col*k/(numchirp*numt));
        end
    end
    
    % range-angle
    N=numty*numry;
    range_angle=squeeze(RDMcom(:,10,:));
    azi_intensity=20*log10(abs(fftshift(fft(range_angle(:,1:N),N,2),2)));
    ele_intensity=20*log10(abs(fftshift(fft(range_angle(:,N+1:2*N),N,2),2)));
    figure
    imagesc(asind(sin_angle_axis),range_axis,azi_intensity); colormap('jet'); axis xy
    %xlim([-2,2]); 
    colorbar
    ylim=[0 10];
    colormap; %xlim([1 9])
    clim = get(gca,'CLim');
    set(gca, 'CLim', clim(2)+[-20,0]);
    xlabel('Azimuth(degree)','FontSize',16);
    ylabel('Range(m)','FontSize',16)
    set(gca, 'FontSize',16)
    figure
    imagesc(asind(sin_angle_axis),range_axis,ele_intensity); colormap('jet'); axis xy
    %xlim([-2,2])
    colorbar
    colormap; %xlim([1 9])
    clim = get(gca,'CLim');
    set(gca, 'CLim', clim(2)+[-20,0]);
    xlabel('Elevation(degree)','FontSize',16);
    ylabel('Range(m)','FontSize',16)
    set(gca, 'FontSize',16)
    I1=[];
    I2=[];
    for i=(1:size(index,2))
        cloudpoints.range(i)=range_axis(det_rangeindex(i));   % get range
        cloudpoints.velocity(i)=doppler_axis(det_veloindex(i)); % get velocity
        data_angle=squeeze(RDMcom(index(1,i),index(2,i),:));
        angle_azi=abs(fftshift(fft(data_angle(1:N),N,1),1));
        [M,I1(i)]=max(angle_azi);
        cloudpoints.azi_angle(i)=asin(sin_angle_axis(I1(i)));  % get azimuth angle
        angle_ele=abs(fftshift(fft(data_angle(N+1:2*N),N,1),1));
        [M,I2(i)]=max(angle_ele);
        cloudpoints.ele_angle(i)=asin(sin_angle_axis(I2(i)));  % get elevation angle
        cloudpoints.x(i)=cloudpoints.range(i)*cos(cloudpoints.ele_angle(i))*cos(cloudpoints.azi_angle(i))-4;
        cloudpoints.y(i)=cloudpoints.range(i)*cos(cloudpoints.ele_angle(i))*sin(cloudpoints.azi_angle(i))+1;
        cloudpoints.z(i)=cloudpoints.range(i)*sin(cloudpoints.ele_angle(i))+1; 
    end
end


