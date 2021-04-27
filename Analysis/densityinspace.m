clear;
clc;
nodes = 100;
dx1 = 40.77998; % Angstrom
density0 = 0.0590; % 1/A^3
 %%
pos = zeros(nodes,1);
volume = zeros(nodes,1);
time = 80000;%[2000; 6000; 10000; 16000];

for i=1:nodes
   pos(i) = dx1 * (i-40) / 10;
end
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
                    volume(i) = 1.0/A(2,1);
                end
                lines = fgetl(rawi);
                A = sscanf(lines,'%f');
            end
            break;
        end
        lines = fgetl(rawi);
    end
    fclose(rawi);
%volume(69) = volume(68);
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
hold all
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
plot(pos,volume/density0,'LineWidth',3);
%line([0 200],[density0 density0],'Color','red','LineStyle','--');
%ylabel('Density (1/A^3)','FontWeight','bold','Fontsize',22);
ylabel('Normalized Density','FontWeight','bold','Fontsize',22);
xlabel('Position (nm)','FontWeight','bold','Fontsize',22);
set(gca,'FontSize',20);
%xlim([0, 200]);
%ylim([0.05, 0.065]);
%xlim([-22, 210]);
%ylim([15 25]);
%legend('Lattice','Electron');
%title('OTM step, 0.3 J/cm^2');
%hold off
