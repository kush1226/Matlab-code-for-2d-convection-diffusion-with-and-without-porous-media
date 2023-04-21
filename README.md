# Matlab-code-for-2d-convection-diffusion-with-and-without-porous-media
MATLAB CODE FOR 2D CONVECTION-DIFFUSION
L = 5; % length of rod
n_points = 101; % n_points = nx = ny
h = L/(n_points-1); % h = dx = dy = L/(n_points-1)
x = 0:h:L;
y = 0:h:L;

T(1, 1:n_points) = 0;
T(1:n_points, 1) = 25;
T_new(1, 1:n_points) = 100;
T_new(1:n_points, 1) = 0;

rho = 1; %density of material
u = 1; %velocity in x direction
v = 1; %velocity in y direction
gamma = 0.05; 

error = 1;
itr = 0;
while error > 1e-7
    for i = 2:n_points-1
	    for j = 2:n_points-1
		    aE = gamma/h - rho*u*h/2;
		    aW = gamma/h + rho*u*h/2;
		    aN = gamma/h - rho*v*h/2;
		    aS = gamma/h + rho*v*h/2;
		    aP = aE + aW + aS + aN + 0;
		    T_new(i,j) = (aE*T(i+1,j) + aW*T(i-1,j) + aN*T(i, j-1) + aS*T(i, j+1))/aP;
	    end
    end
    itr = itr + 1;
    error = 0;
    for i = 2:n_points-1
        for j = 2:n_points-1
            error = error + abs(T(i,j)-T_new(i,j));
        end
    end
    T = T_new;
end


x_L = ((1:n_points)-1).*h;
y_L = 1-((1:n_points)-1).*h;
[X,Y] = meshgrid(x_L, y_L);
contourf(X,Y,T,20);
colorbar;

figure;
plot(1-y,T(:, (n_points+1)/2),'--o');
 









MATLAB CODE FOR 2D CONVECTION-DIFFUSION WITH POROUS MEDIA
L = 5; % length of rod
n_points = 101; % n_points = nx = ny
h = L/(n_points-1); % h = dx = dy = L/(n_points-1)
x = 0:h:L;
y = 0:h:L;

T(1, 1:n_points) = 0;
T(1:n_points, 1) = 25;
T_new(1, 1:n_points) = 100;
T_new(1:n_points, 1) = 0;

rho = 1; %density of material
u = 1; %velocity in x direction
v = 1; %velocity in y direction
mu=1e-2; %dynamic viscosity
k=1; %permiability of porous media
gamma = 0.05;

error = 1;
itr = 0;
while error > 1e-7
    for i = 2:n_points-1
	    for j = 2:n_points-1
		    aE = gamma/h - rho*u*h/2;
		    aW = gamma/h + rho*u*h/2;
		    aN = gamma/h - rho*v*h/2;
		    aS = gamma/h + rho*v*h/2;
		    aP = aE + aW + aS + aN + 0;
		    T_new(i,j) = (aE*T(i+1,j) + aW*T(i-1,j) + aN*T(i, j-1) + aS*T(i, j+1))/aP+(mu*u)*h*h/k+(mu*v)*h*h/k;
	    end
    end
    itr = itr + 1;
    error = 0;
    for i = 2:n_points-1
        for j = 2:n_points-1
            error = error + abs(T(i,j)-T_new(i,j));
        end
    end
    T = T_new;
end


x_L = ((1:n_points)-1).*h;
y_L = 1-((1:n_points)-1).*h;
[X,Y] = meshgrid(x_L, y_L);
contourf(X,Y,T,20);
colorbar;

figure;
plot(1-y,T(:, (n_points+1)/2),'--o');
 


