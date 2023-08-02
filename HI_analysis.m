clear all
close all
%% Read ground
test = csvread('ground.csv',0,6);
test = test(1,:);
d = (10.^(test./10))./1000;
n = length(test);
x = linspace(1419.5,1421.5,n);
%v = 299792.458*(1 - (x./1420.4057511));
%  plot(x,d,'LineWidth',2)
% grid on
Pg = d;
clear d n test;
%% Read sky
test = csvread('sky.csv',0,6);
test = test(1,:);
d = (10.^(test./10))./1000;
n = length(test);
x = linspace(1419.5,1421.5,n);
%v = 299792.458*(1 - (x./1420.4057511));
% plot(x,d,'LineWidth',2)
%grid on
Psky = d;
clear d n test;

%% Read source
test = csvread('seventy.csv',0,6);
test = test(1,:);
d = (10.^(test./10))./1000;
n = length(test);
x = linspace(1419.5,1421.5,n);
%v = 299792.458*(1 - (x./1420.4057511));
%  plot(x,d,'LineWidth',2)
% grid on
Psource = d;
clear d n test;

%%
figure

plot(x,Pg,'LineWidth',2)
hold on
grid on
plot(x,Psky,'LineWidth',2)
plot(x,Psource,'LineWidth',2)

xlabel('Frequency (MHz)')
ylabel('Power (W)')
legend ('Pg','Psky','Psource')
set(gca,'FontSize',15)
%%
Tsky = 5;
%Tsky = 5*(1/cosd(90-65));
Tg = 300;


%% Get Radiometer temperature and source temperature
y = Pg./Psky;
Tr = (y*Tsky - Tg)./(1 - y);
Tr1 = Tr;
x_new = 1:1:length(Tr);
figure
plot(x_new,Tr,'LineWidth',2)
grid on
%%
% enter approximate x y coordinates for two points excluding the HI line
x1 = 192;
y1 = 83.7181;
x2 = 245;
y2 = 83.4529;

% Compute the best-fit line
coefficients = polyfit([x1, x2], [y1, y2], 1);
p1 = coefficients (1);
p2 = coefficients (2);

yfit = p2+x_new.*p1;

% replacing mod h1 line portion with fit line
Tr(:,x1:x2) = yfit(:,x1:x2); % new Tr1
% figure
% plot(v,Tr,'LineWidth',2)
% legend ('Tr')

%% get Tsource and Tbrightness
y1 = Psource./Pg;
Tsource = y1.*(Tg + Tr) - Tr;
Tb = Tsource - mean(Tsource(:,1:3));


%%
figure

plot(x,Tr1,'LineWidth',2)
hold on
grid on
plot(x,Tr,'LineWidth',2)
plot(x,Tsource,'LineWidth',2)
plot(x,Tb,'LineWidth',2)
legend T_{r} T_{fit} T_{source} T_{brightness}
xlabel('Frequency (MHz)')
ylabel('Temperature (K)')
hold off
set(gca,'FontSize',15)

%% Velocity correction
% get v_corr from python code
v = 299792.458*(1 - (x./1420.4057511));
v_corr = -3.3602; % from python script
v_new = v + v_corr; % LSR corrected velocity
figure

plot(v_new,Tb,'LineWidth',2)
legend T_{brightness}
xlabel('V_{LSR} (km/s)')
ylabel('Temperature (K)')
grid on
set(gca,'FontSize',15)

%
%% open curve fitting tool fitting
% v and Tb are the x and y data
% save the sfit file . Here the file is saved as fitting.sfit. The code in
% following section will not work without the .sfit file.
curveFitter

%% all gaussians
filename = "fitting.sfit"; % loading the .sfit file. put the name of the saved .sfit file here.

vd1 = load(filename, '-mat');

parameters = vd1.savedSession.AllFitdevsAndConfigs{1, 1}.Fitdev.Output.numparam;
l = (parameters)/3;
a = zeros(l,1);
b = zeros(l,1);
c = zeros(l,1);
coeff = coeffvalues(vd1.savedSession.AllFitdevsAndConfigs{1, 1}.Fitdev.Fit  );
k = 1;
for j = 1:3:parameters

    a(k) = coeff(j); % peak amplitude
    b(k) = coeff(j+1); % peak velocity
    c(k) = coeff(j+2); % peak width
    k = k + 1;
end
clear j k

f = [];
%k = 1;
residuals = vd1.savedSession.AllFitdevsAndConfigs{1, 1}.Fitdev.Output.residuals;
Legend=cell(length(a),1);
figure
subplot(2,1,1)
for j = 1:1:length(a)
    k = a(j).*exp(-((v_new-b(j))/c(j)).^2);
    f = [f;k];
    %k = k+1;
    plot(v_new,k,'LineWidth',2)
    hold on
    Legend{j}=strcat('Gauss', num2str(j));
    
end

plot(v_new,Tb,'LineWidth',3)

f_fit = sum(f,1);

plot(v_new,f_fit,'LineWidth',2)
Legend{end+1} = 'Tb';
Legend{end+1} = 'fitted model';
legend(Legend)
xlabel('V_{LSR} (km/s)')
ylabel('Temperature (K)')
grid on
set(gca,'FontSize',15)

subplot(2,1,2)

plot(v_new,residuals,'LineWidth',2,'Color','black')
xlabel('V_{LSR} (km/s)')
ylabel('T_{resudial} (K)')
legend('Residuals')
grid on
set(gca,'FontSize',15)

%% max velocity
Tb_sigma = vd1.savedSession.AllFitdevsAndConfigs{1, 1}.Fitdev.Goodness.rmse;
sigma_value = 3;

[xs,ind] = sort(v_new,'desc');
ys = f_fit(ind);

sigma_value = 3;

ind_ys = find(ys>=(sigma_value*Tb_sigma));
ind_sel = ind(ind_ys(1));
vel = v_new(ind_sel);
y_vel = f_fit(ind_sel);
v2 = vel;

ind_ys1 = find(ys>=((sigma_value+1)*Tb_sigma));
ind_sel1 = ind(ind_ys1(1));
vel1 = v_new(ind_sel1);
y_vel1 = f_fit(ind_sel1);
v_upper = vel1;

ind_ys2 = find(ys>=((sigma_value-1)*Tb_sigma));
ind_sel2 = ind(ind_ys2(1));
vel2 = v_new(ind_sel2);
y_vel2 = f_fit(ind_sel2);
v_lower = vel2;


figure
plot(v_new,f_fit,'b-','linewidth',3)
hold on
scatter(vel,y_vel,50,'ro','filled')
scatter(vel1,y_vel1,50,'bo','filled')
scatter(vel2,y_vel2,50,'ko','filled')
box on
set(gca,'fontsize',20)
xlabel('V_{LSR} (Km/s)')
ylabel('Brightness Temperature (K)')
grid on
legend ('Fitted model','3 sigma value','3+1 sigma value','3-1 sigma value')
