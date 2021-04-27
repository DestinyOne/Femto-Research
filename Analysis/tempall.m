clear;
dx1 = 40.77998;
K2eV = 8.621738e-5; 
%%
posdim = 100;
timedim = 50;
pos = zeros(posdim,1);
timee = 1:1:timedim;
Data = zeros(timedim, posdim) + 300;
time = timedim*1000;

for i=1:posdim
   pos(i) = dx1 * (i-40) / 10;
end
posmax = pos(1);
posmin = pos(posdim);
target = 'Ta';
raw = [target 'out.txt'];
rawi = fopen(raw,'r');

lines = fgetl(rawi);
while ischar(lines)
    A = sscanf(lines,'%f');

    if isequal(size(A),[1 1]) && A(1,1) <= time
      ti = A(1,1)/1000;
      % disp(ti);
      fgetl(rawi);
      lines = fgetl(rawi);
      A = sscanf(lines,'%f');
      while isequal(size(A),[2 1])
        i = A(1,1);
        if i >= 1 && i <= posdim    
          posmin = min(posmin,pos(i));
          posmax = max(posmax,pos(i));
          Data(ti,i) = A(2,1);
        end
        lines = fgetl(rawi);
        A = sscanf(lines,'%f');
      end
      % break;
    end
    lines = fgetl(rawi);
end
fclose(rawi);

%% Plot results
% set(gca,'FontSize',30,'FontWeight','bold');
% ylabel('Temperature (K)','FontWeight','bold','Fontsize',30);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.13; b = 0.12; 
c = 0.7; d = 0.8;

figure( 'Position', [200 50 1100 900], 'InvertHardcopy', 'off', 'Color', 'w' ); 
colormap( jet(512) );
han = axes( 'Units', 'normalized', 'Position', [a b c d],...
            'Xlim', [timee(1), timee(timedim)],... 
            'YLim', [-30, 20], 'Ydir', 'reverse',...
            'FontName', 'Times', 'FontSize', 24, 'FontWeight', 'bold', 'Layer', 'top','LineWidth',1.5);

grid on;
grid minor;
hold on;
xlabel( 'Time (ps)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 32);
ylabel( 'Position (nm)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 32); 
text( 0.10, 0.90, sprintf( '%5.0f ps', 10), 'Units', 'normalized',  'HorizontalAlignment', 'center', 'FontName', 'Times', 'FontSize', 50, 'FontWeight', 'bold', 'Color', 'r' );

set(colorbar, 'Position', [0.89 0.12 0.02 0.78], 'FontName', 'Times',...
  'FontSize', 24, 'FontWeight', 'bold' );
caxis([300 9500])
% 'Dirction', 'reverse'
hold on;

pcolor(timee, pos, Data');
shading interp;
hold on;  

text( 1.04, 0.5, 'Temperature (K)', 'Units', 'normalized', 'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'FontName', 'Times', 'FontSize', 32, 'FontWeight', 'bold', 'Color', 'k' );

return

