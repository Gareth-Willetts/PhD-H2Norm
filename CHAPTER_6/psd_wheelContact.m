V = 55; Av = 2.5*10^(-7);
S = 2*pi*V*Av;

x = 0:0.1:200;
y = zeros(1,length(x));
y(x < 40*pi) = S;
plot(x,y, 'LineWidth', 2); grid on;
xlabel("Frequency (rad/s)"); ylabel({'Power Spectral Density of Leading'; 'Wheel Vertical Contact Velocity'});
set(gca, 'FontSize', 14)