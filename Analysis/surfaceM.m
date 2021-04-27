%clear;
%clc;
dx = 40.77998;
 %%
mystep = 500;
time = zeros(mystep,1);
pos = zeros(mystep,1);
vel = zeros(mystep,1);

timed = [0 14.58 19.76 29.63 49.86 74.79 99.72 199.88 299.51 500.00];
posd = [0 32.45 52.81 85.90 119.62 132.34 138.07 288.25 400.14 517.16];

timee = [0 50 75 100 200];
pose = [0 55.28 101.04 81.58 143.75];


for i = 1:mystep
  time(i) = i;
end

rawi = fopen('surface.txt','r');

lines = fgetl(rawi);
while ischar(lines)
    A = sscanf(lines,'%f');
    if A(1,1) == 0
      poso = A(2,1)*dx*150;
    end  
    for i = 1:mystep
      if A(1,1) == time(i)
        pos(i) =  poso-A(2,1)*dx*150;
        break;
      end
    end
    lines = fgetl(rawi);
end
fclose(rawi);

vel(1) = pos(1)/10;

for i = 2:mystep
  vel(i) = (pos(i)-pos(i-1))/10;
end

%% Plot results
%figure
%yyaxis left
%subplot(2,2,3);
%plot(pos,Ta,'-','LineWidth',3);
%plot(time,pos/10,'LineWidth',3);
plot(time,pos/10,'LineWidth',3);
%plot(time,vel,'LineWidth',2);
%legendInfo{i} = [int2str(time/1000) ' ps']; 
%legend([int2str(time/1000) ' ps']);%, '6 ps', '10 ps', '16 ps');
%legend(legendInfo);
hold on;
%scatter(timed, posd/10, 'o', 'MarkerFaceColor', 'r');
%scatter(timed, posd/10, 30,'o', 'MarkerFaceColor', 'r');
plot(timed,posd/10,'-o','LineWidth',3);
%plot(timee,pose/10,'-d','LineWidth',3);
%scatter(timee, pose/10, 30,'o', 'MarkerFaceColor', 'b');
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
set(gca,'FontSize',30,'FontWeight','bold');
xlabel('Time (ps)','FontWeight','bold','Fontsize',30);
ylabel('Displacement (nm)','FontWeight','bold','Fontsize',30);
%ylabel('Velocity (km/s)','FontWeight','bold','Fontsize',22);
grid('on');
xlim([0, 500]);
%xlim([-22, 210]);
%xlim([0, 200]);
%legend('Lattice','Electron');
%title('OTM step, 0.3 J/cm^2');
%title('(c)','Fontsize',20);
%hold off;
