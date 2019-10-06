% Cubic Quintic Complex Ginzburg Landau
% Eqations

clear all;
close all;
clc;

n_sq = @(b,beta,x1, a1, e1) (b * (3*a1*imag(e1)+ real(e1) -  2*(a1^2) * real(e1) ) * x1) ...
    /(a1 * (real(beta) * imag(e1) - imag(beta)* real(e1) ));

mu_sq = @(b,beta,x1, a1, e1) (b * (-3*a1*imag(beta) - real(beta) + 2*(a1^2) * real(beta) ) * x1) ...
    /(a1 * (real(beta) * imag(e1) - imag(beta)* real(e1) ));

r_sq = @(x1, a1, y1) x1 / (a1 * imag(y1));

k2_sq = @(a1, a2, x1, y1) ((a2^2) * x1) / (a1 * imag(y1));

k1 = @(v, a1, y1, x1) (-v + sqrt(a1 * imag(y1) * x1))/ imag(y1);

w = @(r, v, k1, y1, x1) (-r*v - 2 * r * k1 * imag(y1) + 2 * x1);

omega1 = @(v, k1, r_sq, y1) (v * k1 - (r_sq * imag(y1))/4 + k1^2*imag(y1));

omega2 = @(v, k2, r_sq, y2) (- 0.25 * r_sq * imag(y2) + k2^2 * imag(y2) - k2*(v + sqrt(r_sq)* real(y2)));

p_a1 = @(delta1)[-4 * real(delta1), 8 * imag(delta1), 3*real(delta1)];
p_a2 = @(delta2, y2)[4 * (-real(y2)*imag(delta2) + (imag(y2) * real(delta2))), ...
    -8 * (imag(y2)*imag(delta2) + (real(y2) * real(delta2))), ...
    3*real(y2)*imag(delta2) - 3 * imag(y2) * real(delta2)];

P = @(a2, y2) 3 * a2 * imag(y2) - real(y2) + 2 * (a2^2) * real(y2);
y1i = @(a2, y2, beta1, beta2, e1, e2, a1) P(a2, y2) * (real(beta1) * imag(e1) - imag(beta1) * real(e1))/ ...
    (-real(beta1) * real(beta2) + real(e1) * real(e2) + 3*a1 * (-imag(beta1) * real(beta2) + imag(e1)*real(e2)) + (a1^2)*(2*real(beta1)*real(beta2) - 2 * real(e1) * real(e2)));
% TODO: implement logic deciphering
% wether to take the root value less
% than 0 or greater than zero
% a1 = routs(p_a1(delta1));
% a2 = roots(p_a2(delta2, y2));


L = @(b,beta,x1, a1, e1, delta) (b^2)/(8*(a1^2)*(real(beta) * imag(e1) - imag(beta) * real(e1))) * ...
    ((6 * a1 * real(delta) * imag(e1) * real(e1) * x1 * (1-2*(a1^2) + real(delta) * (real(e1)^2)*x1 * (1+4 * (a1^4)))) + ...
    (a1 ^ 2) * (2 * real(beta)^2 * imag(e1)^2 - 4 * imag(beta) * real(beta)*imag(e1)*real(e1) + 2 * imag(beta)^2 * real(e1)^2 + ...
    real(delta)* (9 * imag(e1)^2 - 4*real(e1)^2)*x1));

A = @(x, t, n, r, w, k1, omega1, a1, L, b) (n * exp(0.5 * (r*x+w*t) + 1i*(k1 * x - omega1 * t)))/...
    ((1+b*exp(r*x+w*t) + L * exp(2 * (r*x+w*t)))^(0.5 + 1i *  a1));

B = @(x, t, n, r, w, k2, omega2, a2, L, b) (n * exp(0.5 * (r*x+w*t) + 1i*(k2 * x - omega2 * t)))/...
    ((1+b*exp(r*x+w*t) + L * exp(2 * (r*x+w*t)))^(0.5 + 1i * a2));

% Sample Experiment
beta1 = 1.2;
delta1 = 2.4+0.6i;
e1 = 1.6+0.5i;
x1 = 0.01;
y2 = 0.75-1.1i;
beta2 = 1i;
delta2 = 2.75-i;
e2 = -2;
b = 3;
v = 0;
k1v = 5;
% obtain negative roots
a1 = findNegativeValueInArray(roots(p_a1(delta1)));
a2 = findNegativeValueInArray(roots(p_a2(delta2, y2)));
y1_i = y1i(a2, y2, beta1, beta2, e1, e2, a1);

y1 = 0.4 - 0.7i;
nsq = n_sq(b,beta1,x1, a1, e1); 
rsq = r_sq(x1, a1, y1);
musq = mu_sq(b,beta1,x1, a1, e1); 
rsq = r_sq(x1, a1, y1);
k1val = k1(v, a1, y1, x1);
k2sq = k2_sq(a1, a2, x1, y1);
wval = w(sqrt(rsq), v, k1val, y1, x1);
o1 = omega1(v, k1val, rsq, y1);
o2 =  omega2(v, sqrt(k2sq), rsq, y2);
L1 = L(b,beta1,x1, a1, e1, delta1);
x = 1;
t = 1;
Aval = A(x, t, sqrt(nsq), sqrt(rsq), wval, k1val, o1, a1, L1, b);

%f=1.0e6;
%fs=f*100;
%n=3; 
%t=0:1/fs:3/f;
t = 0:1/40:10;
[X,Y] = meshgrid(-200:1:200, -200:1:200);
t = meshgrid(t);
amplitude = [];
amplitude_b = [];


for i = 1:length(t)
   for j = 1:length(X) 
        if t(1, i) == 0 || X(1, i) == 0
           amplitude(i, j) =  A(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, k1val, o1, a1, L1, b) * (1+0.1 *rand(1));
           amplitude_b(i, j) = B(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, sqrt(k2sq), o2, a2, L1, b) * (1+0.1 * rand(1));
        else
            amplitude(i, j) =  A(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, k1val, o1, a1, L1, b);
            amplitude_b(i, j) = B(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, sqrt(k2sq), o2, a2, L1, b);
        end
   end
end


amplitude = (abs(amplitude));
amplitude_b = (abs(amplitude_b));

subplot(2,1,1)
colormap(cool);
s1 = surf(t,Y,amplitude,'FaceAlpha',0.8); 
s1.EdgeColor = 'none';

%title(sprintf('2D wave equation at t = %1.2f, con sigma = %1.2f y gamma = %1.2f',t(j),sigma, gamma),'Fontsize',11);
title("bright-bright solitary waves for eq(1) (Perturbed)");
xlabel('time'); ylabel('x'); zlabel("|A|^2");
%zlabel(sprintf('u(x,y,t = %1.2f)',t(j)),'Fontsize',11);
%axis ([0 1 0 1 -1 1]);

subplot(2,1,2)
colormap(cool);
s2 = surf(t,Y,amplitude_b,'FaceAlpha',0.8); 
s2.EdgeColor = 'none';

%title(sprintf('2D wave equation at t = %1.2f, con sigma = %1.2f y gamma = %1.2f',t(j),sigma, gamma),'Fontsize',11);
title("bright-bright solitary waves for eq(2) (Perturbed)");
xlabel('time'); ylabel('x'); zlabel("|A|^2");


%%
for i = 1: length(t)
    if t(1,i) == 0 || X(1,i) == 0
        rng('shuffle');
        amplitude = [amplitude,  A(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, k1val, o1, a1, L1, b) * (1+0.7 *rand(1))];
        amplitude_b = [amplitude_b,  B(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, sqrt(k2sq), o2, a2, L1, b) * (1+0.1 * rand(1))];
    else
        amplitude = [amplitude,  A(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, k1val, o1, a1, L1, b)];
        amplitude_b = [amplitude_b,  B(X(1,i), t(1,i), sqrt(nsq), sqrt(rsq), wval, sqrt(k2sq), o2, a2, L1, b)];
    end
end
amplitude = meshgrid(abs(amplitude));
amplitude_b = meshgrid(abs(amplitude_b));

figure(1);

subplot(2,1,1)
colormap(cool);
s1 = surf(t,Y,amplitude,'FaceAlpha',0.8); 
s1.EdgeColor = 'none';

%title(sprintf('2D wave equation at t = %1.2f, con sigma = %1.2f y gamma = %1.2f',t(j),sigma, gamma),'Fontsize',11);
title("bright-bright solitary waves for eq(1)");
xlabel('time'); ylabel('x'); zlabel("|A|^2");
%zlabel(sprintf('u(x,y,t = %1.2f)',t(j)),'Fontsize',11);
%axis ([0 1 0 1 -1 1]);

subplot(2,1,2)
colormap(cool);
s2 = surf(t,Y,amplitude_b,'FaceAlpha',0.8); 
s2.EdgeColor = 'none';

%title(sprintf('2D wave equation at t = %1.2f, con sigma = %1.2f y gamma = %1.2f',t(j),sigma, gamma),'Fontsize',11);
title("bright-bright solitary waves for eq(2)");
xlabel('time'); ylabel('x'); zlabel("|A|^2");

function val = findNegativeValueInArray(arr)
    [row, col] = find(arr < 0);
    val = arr(row, col);
end


function val = findPositiveValueInArray(arr)
    [row, col] = find(arr > 0);
    val = arr(row, col);
end
