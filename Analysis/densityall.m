clear;
dx1 = 40.77998;
K2eV = 8.621738e-5; 
%%
posdim = 100;
timedim = 300;
pos = zeros(posdim,1);
timee = 1:1:timedim;
Data = zeros(timedim, posdim);
time = timedim*1000;
density0 = 0.0590; % 1/A^3

for i=1:posdim
   pos(i) = dx1 * (i-40) / 10;
end
posmax = pos(1);
posmin = pos(posdim);
raw = "volume.txt";
rawi = fopen(raw,'r');

lines = fgetl(rawi);
while ischar(lines)
    A = sscanf(lines,'%f');

    if isequal(size(A),[1 1]) && A(1,1) <= time 
      ti = A(1,1)/1000;
      if ti == 0
        ti = 1;
      end
      % disp(ti);
      fgetl(rawi);
      lines = fgetl(rawi);
      
      A = sscanf(lines,'%f');
      while isequal(size(A),[2 1])
        i = A(1,1);
        if i >= 1 && i <= posdim    
          posmin = min(posmin,pos(i));
          posmax = max(posmax,pos(i));
          Data(ti,i) =  1.0/A(2,1);
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = 0.13; b = 0.12; 
c = 0.7; d = 0.8;

figure( 'Position', [200 50 1100 900], 'InvertHardcopy', 'off', 'Color', 'w' ); 
colormap( jet(512) );

han = axes( 'Units', 'normalized', 'Position', [a b c d],...
            'Xlim', [timee(1), timee(timedim)],... 
            'YLim', [-50, 50], 'Ydir', 'reverse',...
            'FontName', 'Times', 'FontSize', 32, 'FontWeight', 'bold', 'Layer', 'top','LineWidth',1.5);

grid on;
grid minor;
hold on;
xlabel( 'Time (ps)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 32);
ylabel( 'Position (nm)', 'FontName', 'Times', 'FontWeight', 'bold', 'FontSize', 32); 
% text( 0.16, 0.90, sprintf( '%5.0f  mJ/cm^2', 100 ), 'Units', 'normalized',  'HorizontalAlignment', 'center', 'FontName', 'Times', 'FontSize', 34, 'FontWeight', 'bold', 'Color', 'r' );
text( 0.13, 0.93, sprintf( '%5.0f  ps', 10 ), 'Units', 'normalized',  'HorizontalAlignment', 'center', 'FontName', 'Times', 'FontSize', 50, 'FontWeight', 'bold', 'Color', 'r' );

colorbar('Position', [0.89 0.12 0.02 0.8], 'FontName', 'Times','FontSize', 32, 'FontWeight', 'bold');
caxis([0 1.1])

hold on;

pcolor(timee, pos, Data' / density0);
shading interp;
hold on;  

text( 1.04, 0.5, 'Normalized Density', 'Units', 'normalized', 'Rotation', 90, 'HorizontalAlignment', 'center', ...
    'FontName', 'Times', 'FontSize', 32, 'FontWeight', 'bold', 'Color', 'k' );

%%
% S-L interface
raw = 'liquid.txt';
rawi = fopen(raw,'r');
posls = zeros(timedim,1);
for i = 1:timedim+1
  lines = fgetl(rawi);
  A = sscanf(lines,'%f');
  ti = A(1,1)/1000;
  if ti == 0
    continue;
  end
  posls(ti) = dx1 * (A(2,1)-40) / 10;
end
fclose(rawi);
plot(timee,posls,'-k','LineWidth',3);

%v = [1.0, 1.0];
%[C, hcNe] = contour(timee, pos, Data' ./ density0, v, 'r','LineWidth', 3.,'LineStyle','--');
%clabel(C, hcNe, 'FontName', 'Times', 'FontSize', 26, 'FontWeight', 'bold','Color','r');

%print ('-dbmp', 'Fig1.bmp');
%close( gcf );
return

