clear;
clc;
%%
multi = 1;
dx1 = 40.77998; % Angstrom
density0 = 0.0590; % 1/A^3
dV = 1*1*1*dx1*dx1*dx1;
nodes = 200;
% nodes = nodes + 1;
raw = 'liquid.txt';
rawi = fopen(raw,'r');
pos = zeros(nodes,1);
num = zeros(nodes,1);
thick = zeros(nodes,1);
time = zeros(nodes,1);
linea = fgetl(rawi);
count = 1;
while ischar(linea)   
    A = sscanf(linea,'%f');
    time(count) = A(1,1)/1000;
    pos(count) = dx1 * (A(2,1)/multi-40) / 10;
    num(count) = A(3,1);
    thick(count) = num(count)/density0/dx1/dx1/3/3/10; %nm

    linea = fgetl(rawi);
    count = count + 1;
end
fclose(rawi);

%% Plot results
% figure
%subplot(1,2,2);
plot(time,thick,'-r','LineWidth',3);
set(gca,'FontSize',30,'FontWeight','bold');
ylabel('Melted thickness (nm)','FontWeight','bold','Fontsize',30);
xlabel('Time (ps)','FontWeight','bold','Fontsize',30);
xlim([2 nodes]);

grid('on');
grid minor;
hold on;
%subplot(1,2,1);
%plot(time,pos,'LineWidth',3);
%plot(time,Te,'LineWidth',3);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
%set(gca,'FontSize',30,'FontWeight','bold');
%ylabel('S-L interface position (nm)','FontWeight','bold','Fontsize',30);
%xlabel('Time (ps)','FontWeight','bold','Fontsize',30);
%grid('on');
%grid minor;
%xlim([5 300]);

%grid('on');
%legend('Lattice','Electron','Location','Southeast');
%title('(c)','Fontsize',20);
%title('2.122\mu Au (Modified Ke)','Fontsize',28);
%hold off
