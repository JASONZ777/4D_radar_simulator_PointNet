% figure
% t=[1:8:48;1:8:48]; %x,y
% r=[1:1:8;1:1:8];
% plot(t(1,:),zeros(1,6),'xb', 'MarkerSize',8,'LineWidth',1)
% hold on
% plot(r(1,:),zeros(1,8)+1,'+r',  'MarkerSize',8,'LineWidth',1)
% hold on
% vy=zeros(6,8);
% vz=zeros(6,8);
% %ty
% for k=1:6
%     vy(k,:)=(t(1,k)+r(1,:))./2;
%     vz(k,:)=(1+r(2,:))./2;
%     plot(zeros(1,8)+(t(1,k)+1)./2,vz(k,:),'.c','MarkerSize',8,'LineWidth',1)
%     hold on
% end
% vy=reshape(vy,1,48);
% plot(vy,zeros(1,48)+0.5,'.c','MarkerSize',8,'LineWidth',1)
% plot(zeros(1,6),t(2,:),'xb', 'MarkerSize',8,'LineWidth',1)
% hold on
% plot(zeros(1,8)+1,r(2,:),'+r', 'MarkerSize',8,'LineWidth',1)
% hold on
% vy=zeros(6,8);
% vz=zeros(6,8);
% %tz
% for k=1:6
%     vy(k,:)=(1+r(1,:))./2;
%     vz(k,:)=(t(2,k)+r(2,:))./2;
%     plot(vy(k,:),zeros(1,8)+(t(2,k)+1)./2,'.c','MarkerSize',8,'LineWidth',1)
% end
% vz=reshape(vz,1,48);
% plot(zeros(1,48)+0.5,vz,'.c','MarkerSize',8,'LineWidth',1)
% xlabel("Distance(x wavelengths)")
% ylabel("Distance(x wavelengths)")
% legend("Tx","Rx","Virtual array")
tx=[1,5];
rx=[1,2,3,4];
vx=[1:0.5:4.5];
plot(tx,zeros(1,2),'xb', 'MarkerSize',12,'LineWidth',2)
hold on
plot(rx,zeros(1,4)+2,'+r',  'MarkerSize',12,'LineWidth',2)
hold on
plot(vx,zeros(1,8)+1,'.c','MarkerSize',12,'LineWidth',2)
xlabel("Distance(x wavelengths)")
ylabel("Distance(x wavelengths)")
legend("Tx","Rx","Virtual array")

