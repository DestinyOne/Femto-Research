clear;
clc;
%%
raw = 'Taout.txt';
rawi = fopen(raw,'r');
mystep = 300;
Ta = zeros(mystep,1);
time = zeros(mystep,1);
i = 1;
line = fgetl(rawi);
while ischar(line)   
    A = sscanf(line,'%f');

    if isequal(size(A),[1 1])
        time(i) = A(1,1)/1000;
        fgetl(rawi);
        line = fgetl(rawi);
        A = sscanf(line,'%f');
        Ta(i) = A(2,1);
        i = i+1;
        if i > mystep
          break;
        end
    end
    
    line = fgetl(rawi);
end
fclose(rawi);
%%
raw = 'Teout.txt';
rawi = fopen(raw,'r');
Te = zeros(mystep,1);
i = 1;
line = fgetl(rawi);
while ischar(line)   
    A = sscanf(line,'%f');

    if isequal(size(A),[1 1])
        time(i) = A(1,1)/1000;
        fgetl(rawi);
        line = fgetl(rawi);
        A = sscanf(line,'%f');
        Te(i) = A(2,1);
        i = i+1;
        if i > mystep
          break;
        end
    end
    
    line = fgetl(rawi);
end
fclose(rawi);
%% Plot results
% figure
%subplot(2,2,3);
plot(time,Ta./1000,'LineWidth',2);
hold on;
%plot(time,Te./1000,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
ylabel('Temperature (kK)','FontWeight','bold','Fontsize',22);
xlabel('Time (ps)','FontWeight','bold','Fontsize',22);
%ylim([0 50]);
grid('on');
set(gca,'FontSize',20);
legend('Lattice','Electron','Location','north');
%title('2.122\mu Au (Modified Ke)','Fontsize',28);
title('(c)','Fontsize',20);
%hold off;
