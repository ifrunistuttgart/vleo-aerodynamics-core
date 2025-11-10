
clear;
global v centroid areas

v = [1;0;0];
n = [0;0;-1];
centroid = [0;0;0];
areas = [1];
delta = 0;
[F,T] = newModel(areas,n,centroid,v,delta,1,"C:\Users\Jan_L\OneDrive\Dokumente\Studium\Bacherlorarbeit\vleo-aerodynamics-tool\cl_cd_cVAE_A01_flat_and_bird.csv");

% Starting point for both vectors
x = 0; y = 0; z = 0;

% Plot the vectors
figure
hold on
grid on
axis equal
xlabel('X')
ylabel('Y')
zlabel('Z')
view(3)

% Plot vector v
quiver3(x,y,z, v(1), v(2), v(3), 0, 'r', 'LineWidth', 2)

% Plot vector n
normal = quiver3(x, y, z, n(1), n(2), n(3), 0, 'b', 'LineWidth', 2);

% force vector
f = quiver3(x,y,z,1,1,1,0,"g","LineWidth",2);
plot_aero_force(f,n,[0]);

%AOA display
txt = uicontrol("Style","Text", ...
    "String","AOA: 0", ...
    "Units", "normalized", ...
    "Position",[0.7, 0.4, 0.2, 0.05]);
%slider
sld = uicontrol("Style","slider", ...
  "Min",0,"Max",180,"Value",0, ...
  "Callback",@(src,event) slider_callback(src,txt,normal,n,f));

function slider_callback(slider,txt,normal,n,f)
    global v
    val = get(slider,"value");
    set(txt,"String",sprintf("AOA: %.2f",val))
    R = [cosd(val) 0 sind(val);0 1 0; -sind(val) 0 cosd(val)];
    n_new = R*n;
    delta__rad = acos(dot(-v, n_new)/(norm(v)*norm(n_new)));
    set(normal,"udata",n_new(1),"vdata",n_new(2),"wdata",n_new(3))
    plot_aero_force(f,n_new,[delta__rad])
end

function plot_aero_force(f,n,delta)

    global areas v centroid 
    [F,T] = newModel(areas,n,centroid,v,delta,1,"C:\Users\Jan_L\OneDrive\Dokumente\Studium\Bacherlorarbeit\vleo-aerodynamics-tool\cl_cd_cVAE_A01_flat_and_bird.csv");
    set(f,"udata",F(1),"vdata",F(2),"wdata",F(3))
end


legend('realtive velocity', 'surfcae normal',"aero force")
title('3D Vectors v and n from Origin')