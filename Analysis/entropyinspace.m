clear;
clc;
multi = 1;
nodes = 100*multi;
dx1 = 40.77998; % Angstrom
density0 = 0.0590; % 1/A^3
dV = 1*1*1*dx1*dx1*dx1;
%%
pos = zeros(nodes,1);
entropy = zeros(nodes,1);
density = zeros(nodes,1);
Ta = zeros(nodes,1);
ds = zeros(nodes,1);
time = 16000;%[2000; 6000; 10000; 16000];

for i=1:nodes
   pos(i) = dx1 * (i/multi-40) / 10;
   %pos(i) = i;
end
raw = "entropy.txt";
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
            if i >= 1 && i <= nodes                    
                entropy(i) = A(2,1);
            end
            line = fgetl(rawi);
            A = sscanf(line,'%f');
        end
        break;
    end
    line = fgetl(rawi);
end
fclose(rawi);
for i=2:nodes
    ds(i) = (entropy(i) - entropy(i-1))/dx1*multi; 
    edge = 18;
    if i <= 2 + edge || i >= nodes -edge % Eliminate edge effects
        ds(i) = 0;
    end    
end
     
%%
raw = "volume.txt";
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
                density(i) = 1.0/A(2,1);
            end
            lines = fgetl(rawi);
            A = sscanf(lines,'%f');
        end
        break;
    end
    lines = fgetl(rawi);
end
fclose(rawi);
density(69) = density(68);
density(70) = density(68);
if time == 0
    density(20) = density0;
end
%%
raw = "Taout.txt";
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
            if i >= 1 && i <= 80                    
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
%if strcmp(target,'Ta' ) 
%    plot(pos,density,'LineWidth',3);
%else
%    plot(pos,density,'--','LineWidth',3);
%end
%right = 600+time/1000*2;
%legendInfo{i} = [int2str(time/1000) ' ps']; 
%legend([int2str(time/1000) ' ps']);%, '6 ps', '10 ps', '16 ps');
%legend(legendInfo);
%hold all
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
%plot(pos,entropy,'-m','LineWidth',3);
plot(pos,ds,'LineWidth',3);
%plot(pos,density*1000,'-','LineWidth',3);
ylabel('Entropy gradient, ds/dx (K_b/A)','color','black','FontWeight','bold','Fontsize',22);
%yyaxis left;
%plot(pos, Ta,'-k','LineWidth',3);
%ylabel('Lattice Temperature (K)','color','black','FontWeight','bold','Fontsize',22);
%plot(pos, density*1.e3,'-m','LineWidth',3);
%ylabel('Density (nm^-^3)','color','black','FontWeight','bold','Fontsize',22);
%set(gca, 'YScale', 'log', 'FontSize',20);
hold on;
%yyaxis right;
%ylabel('Entropy per atom (k_b)','color','black','FontWeight','bold','Fontsize',22);
%plot(pos, entropy,'--k','LineWidth',3);
%ylabel('Entropy per atom','FontWeight','bold','Fontsize',22);
xlabel('Position (nm)','FontWeight','bold','Fontsize',22);
set(gca,'FontSize',20);
%xlim([-22, 210]);
%xlim([0, 200]);
%legend('Lattice','Electron');
title('1 J/cm^2, 201.56 nm Au (500K)');
%hold off
