clear;
clc;
dx1 = 40.77998;
%%
nodes = 100;
pos = zeros(nodes,1);
Pi = zeros(nodes,1);
time = 1000;%[2000; 6000; 10000; 16000];

for i=1:nodes
    pos(i) = dx1 * (i-40) / 10;
end

%% Ion pressure
raw = 'pressure.txt';
rawi = fopen(raw,'r');

line = fgetl(rawi);
while ischar(line)
    A = sscanf(line,'%f');

    if isequal(size(A),[1 1]) && A(1,1) == time
        fgetl(rawi);
        line = fgetl(rawi);
        A = sscanf(line,'%f');
        while isequal(size(A),[2 1])
            i = A(1,1);                 
            Pi(i) = A(2,1)/10000;
            line = fgetl(rawi);
            A = sscanf(line,'%f');
        end
        break;
    end
    line = fgetl(rawi);
end
fclose(rawi);

%% Electron Pressure

%% Plot results
%figure
%yyaxis right;
plot(pos,Pi,'-','LineWidth',3);
%right = 600+time/1000*2;
%legendInfo{i} = [int2str(time/1000) ' ps']; 
%legend([int2str(time/1000) ' ps']);%, '6 ps', '10 ps', '16 ps');
%legend(legendInfo);
hold all
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
ylabel('Stress (GPa)','FontWeight','bold','Fontsize',22);
xlabel('Position (nm)','FontWeight','bold','Fontsize',22);
set(gca,'FontSize',20);
%xlim([-22, 213]);
%xlim([0, 200]);
%legend('Lattice','Electron');
%title('OTM step, 0.3 J/cm^2');
%hold off
