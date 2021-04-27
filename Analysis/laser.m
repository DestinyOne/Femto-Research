clear;
clc;
dx1 = 40.77998;
echarge = 1.6021765e-19;
%%
pos = zeros(80,1);
length = zeros(80,1);
intense = zeros(80,1);
time = 80;%[2000; 6000; 10000; 16000];

for i=1:80
   pos(i) = dx1 * (i-20) / 10;
end

raw = 'Laserout.txt';
rawi = fopen(raw,'r');

lines = fgetl(rawi);
while ischar(lines)
    A = sscanf(lines,'%f');

    if isequal(size(A),[2 1]) && A(1,1) == time
        fgetl(rawi);
        lines = fgetl(rawi);
        A = sscanf(lines,'%f');
        while isequal(size(A),[3 1])
            i = A(1,1);
            if i >= 1 && i <= 80
                length(i) = A(2,1);
                intense(i) = A(3,1)*echarge*1.e28;%/(1-0.8583);
            end
            lines = fgetl(rawi);
            A = sscanf(lines,'%f');
        end
        break;
    end
    lines = fgetl(rawi);
end
fclose(rawi);
%%
timeme = zeros(319,1);
reflect = zeros(319,1);
rawi = fopen(raw,'r');
lines = fgetl(rawi);
i = 1;
while ischar(lines)
    A = sscanf(lines,'%f');
    if isequal(size(A),[2 1])
        timeme(i) = A(1,1);
        reflect(i) = A(2,1);
        i = i+1;
    end
    lines = fgetl(rawi);
end
fclose(rawi);

all = 0;
fall = 0;

width = 80;
for i = 1:319
  fall = fall + sqrt(4*log(2)/pi)*exp(-4*log(2)*(timeme(i)-2*width)^2/width^2);
  all = all + sqrt(4*log(2)/pi)*exp(-4*log(2)*(timeme(i)-2*width)^2/width^2)*(1-reflect(i));
end

Fluence = 0.6; %J/cm^2
%I0 = 1.e+13; %W/cm^2
dt = 1.e-15; %s
I0 = Fluence/dt/width;
fabs = all * I0*dt % J/cm^2
R = 1 - all/fall;

plot(timeme,reflect,'LineWidth',3);
ylabel('Reflectivity','FontWeight','bold','Fontsize',22);
xlabel('Time (fs)','FontWeight','bold','Fontsize',22);
%xlim([0, 170]);
%% Plot results
% figure
%hold on
%plot(pos,intense,'LineWidth',3);
hold on
%line([1 200],[intense(20)/2.718281828 intense(20)/2.718281828],'Color','red','LineStyle','--');
% axis([x(1) x(lines) 0 (floor(max(y))+2)])
%ylabel('Intensity (W/cm^2)','FontWeight','bold','Fontsize',22);
%xlabel('Position (nm)','FontWeight','bold','Fontsize',22);
%xlim([-0, 60]);
grid('on');
set(gca,'FontSize',20);
%title('OTM step (During laser)','Fontsize',28);
