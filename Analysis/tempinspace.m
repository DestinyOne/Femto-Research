clear;
clc;
dx1 = 40.77998;
%%
nodes = 100;
pos = zeros(nodes,1);
Ta = zeros(nodes,1);
time = 30000;%[2000; 6000; 10000; 16000];

for i=1:nodes
   pos(i) = dx1 * (i-40) / 10;
end
target = 'Ta';
raw = [target 'out.txt'];
rawi = fopen(raw,'r');

lines = fgetl(rawi);
while ischar(lines)
    A = sscanf(lines,'%f');

    if isequal(size(A),[1 1]) && A(1,1) == time
        fgetl(rawi);
        lines = fgetl(rawi);
        A = sscanf(lines,'%f');
        while isequal(size(A),[2 1])
            i = A(1,1);
            if i >= 1 && i <= nodes                  
                Ta(i) = A(2,1);
            end
            lines = fgetl(rawi);
            A = sscanf(lines,'%f');
        end
        break;
    end
    lines = fgetl(rawi);
end
fclose(rawi);

%% Plot results
%figure
%yyaxis left
if strcmp(target,'Ta' ) 
    plot(pos,Ta,'-b','LineWidth',3);
else
    plot(pos,Ta,'-r','LineWidth',3);
end
%right = 600+time/1000*2;
%legendInfo{i} = [int2str(time/1000) ' ps']; 
%legend([int2str(time/1000) ' ps']);%, '6 ps', '10 ps', '16 ps');
%legend(legendInfo);
hold on;
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
set(gca,'FontSize',30,'FontWeight','bold');
ylabel('Temperature (K)','FontWeight','bold','Fontsize',30);
xlabel('Position (nm)','FontWeight','bold','Fontsize',30);

%xlim([-22, 210]);
%xlim([0, 300]);
%legend('Lattice','Electron');
%title('OTM step, 0.3 J/cm^2');
%hold off
