clear all
clc

x = 0.1:1/22:1;
y = ((1 + 0.6*sin(2*pi*x/0.7)) + 0.3*sin(2*pi*x))/2;

figure(1)
hold on; grid on;
plot(x, y,"o-")
y = SBF(x,y); %funkcija
plot(x, y,"o--")
hold off;
legend("Duota funkcija","Gauta funkcija")

function y = SBF(x,d)
w0 = randn(1);
w1 = randn(1);
w2 = randn(1);

c1 = 0.1892;
c2 = 0.9088;
r1 = 0.1948;
r2 = 0.2101;

eta = 0.1;
eta_c = 0.001;
eta_r = 0.001;

for j=1:20000
    for i=1:20
        F1=exp(-((x(i)-c1)^2)/(2*(r1^2)));
        F2=exp(-((x(i)-c2)^2)/(2*(r2^2)));

        y(i) = F1*w1+F2*w2+w0;

        e = d(i)-y(i);

        w1 = w1 + eta * e * F1;
        w2 = w2 + eta * e * F2;
        w0 = w0 + eta * e;
        
        % c1 = c1 + eta_c * e * w1 * F1 * (x(i) - c1) / (r1^2);
        % c2 = c2 + eta_c * e * w2 * F2 * (x(i) - c2) / (r2^2);
        % r1 = r1 + eta_r * e * w1 * F1 * ((x(i) - c1)^2) / (r1^3);
        % r2 = r2 + eta_r * e * w2 * F2 * ((x(i) - c2)^2) / (r2^3);
        
    end
end
end