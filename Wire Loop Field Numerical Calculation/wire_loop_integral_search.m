format long;
close all;

MU_0 = 1.25663706127E-6;
R = 0.06;

% Composite Simpson Resolution
num_intervals = 100;

x_dom = [0 10];
y_dom = [-0.05 0.05];
z_dom = [-0.05 0.05];
nx = 100;
ny = 10;
nz = 10;
x_mesh = (x_dom(2) - x_dom(1)) / nx;
y_mesh = (y_dom(2) - y_dom(1)) / ny;
z_mesh = (z_dom(2) - z_dom(1)) / nz;
B_x = zeros(nx, ny, nz);
B_y = zeros(nx, ny, nz);
B_z = zeros(nx, ny, nz);

f_x = @(x, y, z, phi) (z * cos(phi) - y * sin(phi) + R) / ((x^2 + y^2 + z^2 + R^2 - (2*R)*(y*cos(phi) + z*sin(phi)))^(3/2));
f_y = @(x, y, z, phi) (x * sin(phi)) / ((x^2 + y^2 + z^2 + R^2 - (2*R)*(y*cos(phi) + z*sin(phi)))^(3/2));
f_z = @(x, y, z, phi) (-x * cos(phi)) / ((x^2 + y^2 + z^2 + R^2 - (2*R)*(y*cos(phi) + z*sin(phi)))^(3/2));

for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            f_pos_x = @(phi) f_x(x_dom(1) + x * x_mesh, y_dom(1) + y * y_mesh, z_dom(1) + z*z_mesh, phi);
            f_pos_y = @(phi) f_y(x_dom(1) + x * x_mesh, y_dom(1) + y * y_mesh, z_dom(1) + z*z_mesh, phi);
            f_pos_z = @(phi) f_z(x_dom(1) + x * x_mesh, y_dom(1) + y * y_mesh, z_dom(1) + z*z_mesh, phi);
            B_x(x, y, z) = ((MU_0 * R) / (4 * pi)) * comp_simpson(f_pos_x, num_intervals, 0, 2*pi);
            B_y(x, y, z) = ((MU_0 * R) / (4 * pi)) * comp_simpson(f_pos_y, num_intervals, 0, 2*pi);
            B_z(x, y, z) = ((MU_0 * R) / (4 * pi)) * comp_simpson(f_pos_z, num_intervals, 0, 2*pi);
        end
    end
end
            
B = zeros(nx, ny, nz);
for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            B(x, y, z) = sqrt(B_x(x, y, z)^2 + B_y(x, y, z)^2 + B_z(x, y, z)^2);
        end
    end
end

xx = linspace(x_dom(1), x_dom(2), nx);
yy = linspace(y_dom(1), y_dom(2), ny);
zz = linspace(z_dom(1), z_dom(2), nz);

figure(1)
plot(xx, B(1:nx, 5, 5));
figure(2)
plot(xx, B_y(1:nx, 1, 1));


% B_mat = zeros(ny, nz);
% for y = 1:ny
%     for z = 1:nz
%         B_mat(y, z) = B(26, y, z);
%     end
% end
% 
% B_mat
% 
% surf(yy, zz, B_mat)

writematrix(B_x, "Bx_raw.txt");
writematrix(B_y, "By_raw.txt");
writematrix(B_z, "Bz_raw.txt");