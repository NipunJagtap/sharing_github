% 'rectification' function is defined to calculate the number of the stages required in the rectification section of the column.
% Firstly, Rmin is calculatd from the thermal condition (q) of the feed-in, like vapor fraction and the temperature of the fluid
% Secondly, the number of stages in the rectification section is assumed to be maximum of 100, but only few are included in the set if 'x' value is greater
% than the feed concentration (x_F);
% Thirdly, for every stage X and Y values are calculated separately.
% Later, the temperature at every stage is calculated in correspondence to the operating pressure of column and molar fraction of light key on each stage.
% Lastly, it gives 3 graphs. One for the comparison of X and Y values at each stage. Whereas, the other figure shows temperature variation across stages.

% NOTE:-
% This function only works for the binary component system distillation.
# The feed is assumed to be at bubble point(q = 0)

% Following are the variables required for the function to run-
% ant_coeff = The 2*3 matrix of antoines coefficient, where A, B, C values of light key is in three cells of first row
% x_F = Mol fraction of light key in feed
% x_D = Mol fraction of light key in distillate
% Temp_min = Boiling point of light key at operating pressure (in Deg C)
% Temp_max = Boiling point of heavy key at operating pressure (in Deg C)
% alpha = Average reltive volatility of solution
% R_by_Rm = Ratio of actual Reflux ratio (R) to the min Reflux ratio (Rmin)
% op_press = Operating pressure of the column (in torr)


function Table = rectifycol(ant_coeff, x_F, x_D, Temp_min, Temp_max, alpha, R_by_Rm, op_press)
q = 0;

% Function for Rmin calculation 
f =@(r) ((r*x_F)/(r*(1-x_F))) - alpha*(-x_D+x_F*(r+1))/((r+1)*(1-x_F)+x_D-1);

% Finding Rmin value using Secant method
k = 2;
x(1) = 1;
x(2) = 1.2;
while abs(f(x(k))) > 0.01;
  x(k+1) = x(k) - f(x(k)) * (x(k) - x(k-1))/(f(x(k)) - f(x(k-1)));
  k = k+1;
endwhile
Rmin = x(k);
R = Rmin * R_by_Rm;
N = 100;
d = 0;
% Getting a sets of X & Y
X =[];
Y = [];
y(1) = x_D;
for j = [1:N]
x(j) = y(j)/(y(j) + (1-y(j))* 2.5);
  if x(j) > x_F;
    X = [X x(j)];
    Y = [Y y(j)];
    d = d + 1;
  endif
j = j+1;
y(j) = (x(j-1) * R / (R+1)) + (x_D/(R+1));
endfor
% Calculation of temperature at each stage
T = [Temp_min:0.1:Temp_max];
t = [];
c = 0;
for x = X;
  for i = T;
    if (x * exp(ant_coeff(1,1) - (ant_coeff(1,2)/(i+273-ant_coeff(1,3)))) + (1 - x) * exp(ant_coeff(2,1) - (ant_coeff(2,2)/(i+273-ant_coeff(2,3)))) - op_press) < 0.001;
      c = i;
    endif
  endfor
  t = [t c];
endfor
%final table of results
Table = [1:d;X;Y;t].';
% plot of X vs N
plot(1:d,Table(:,2))
xlabel('Number of stages'); ylabel('Mol fraction'); title('X vs Stage number');
% plot of Y vs N
figure (3)
plot(1:d,Table(:,3))
xlabel('Number of stages'); ylabel('Mol fraction'); title('Y vs Stage number');
% plot of Temperature vs N
figure (3)
plot(1:d,Table(:,4))
xlabel('Number of stages'); ylabel('Temperature (Deg Celcius)');title('Temperature vs Stage number');
end
%For a trial run, a solution of benzene and toluene is considered
%Code runs with 'Table = rectify_col(Mat, 0.4, 0.99, 80, 110, 2.5, 2, 760)'
% where, Mat = antoines coefficient of benzene and toluene ([15.901	2788.6	52.36; 16.014	3096.5	53.67])