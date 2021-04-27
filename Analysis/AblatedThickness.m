clear;
clc;
nodes = 100;
dx1 = 40.77998; % Angstrom
density0 = 0.0590; % 1/A^3
ds = density0 * dx1 * dx1 *3 *3;
dV = ds * dx1;
 %%
pos = zeros(nodes,1);
density = zeros(nodes,1);
time = 200000;

for i=1:nodes
   pos(i) = dx1 * (i-40) / 10;
end
    raw = "density.txt";
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
                    density(i) = A(2,1);
                end
                line = fgetl(rawi);
                A = sscanf(line,'%f');
            end
            break;
        end
        line = fgetl(rawi);
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
hold all
%plot(time,Te,'LineWidth',2);
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
% plot(pos,density,'LineWidth',3);
plot(density/dV,'LineWidth',3);
%line([0 200],[density0 density0],'Color','red','LineStyle','--');
%ylabel('Density (1/A^3)','FontWeight','bold','Fontsize',22);
ylabel('Normalized Density','FontWeight','bold','Fontsize',22);
xlabel('Position (nm)','FontWeight','bold','Fontsize',22);
set(gca,'FontSize',20);
%hold off
unnum = 0;
start = 43;
for i = start:nodes
    unnum = unnum + density(i);
end
thickness = 204.0413 - unnum/ds/10; % in nm
disp(thickness)


