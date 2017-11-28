import Simplex

close all
clear
clc

n_q = 10; %number of random charges in a 1x1x1 box
Q = rand(n_q, 3) - 0.5; %generate 10 random charges
q = Q(end, :); %unknown charge
Q = Q(1:end-1, :); %known charge
M = [];
%mapping the measurers box
for i = -0.5:0.5:0.5
    for j = -0.5:0.2:0.5
        %X planes
        M = [M; [-0.5, i, j]];
        M = [M; [0.5, i, j]];

        %Y planes
        M = [M; [i, -0.5, j]];
        M = [M; [i, 0.5, j]];

        %Z planes
        M = [M; [i, j, -0.5]];
        M = [M; [i, j, 0.5]];
    end
end
global n_q
global Q
global q
global M

disp('known charge');
disp(q);


%initializing Simplex class
s = Simplex(@cost, {@bound1}, .25, 1e-20, 250);
s.dt = 0; %animation delta time between frames (0 = off)
s.field = 2; %figure subspace of view
s.slices = 10; %15 planes to draw the isolevel maps
s.color = 'green'; %simplex politope color
s.plot = true; %enable the polotting

[value, coordinates, flips, halvings, area] =  s.compute() %compute the algorithm

error = abs(cost(coordinates(1), coordinates(2), coordinates(3)) - cost(q(1), q(2), q(3)));
disp('error');
error

function f = bound1(x, y, z) %sphere bound
    f = -((x+1.5).^2 + (y+1.5).^2 + (z-.8).^2 - 2^2);
end

function f = bound2(x, y, z) %sphere bound
    f = -((x+1.5).^2 + (y+1.5).^2 + (z+.5).^2 - 2^2);
end

function f = get_measure(Q, m) %get the potential from the measure in a certain position
    m = repmat(m, [length(Q(1:end, :)), 1]); %adapt the matrix
    D = sqrt(sum((Q-m).^2, 2)); %norm
    
    f = 1/sum(D); %potential
    f = f/sqrt(3); %normalization factor: 1x1x1 box diagonal 
end

function f = cost(x, y, z)
    global n_q
    global Q
    global q
    global M
    
    %get the measure for each measurer for all the charges
    P = zeros(1, length(M(:, 1))); 
    for j = 1:length(M(:, 1))
        P(1, j) = get_measure([Q; q], M(j, :));  %computes for the known charges plus the unknown charge
    end
    P = sum(P); %sum all the measures
    
    %get the measure for each measurer for all the chargex except the unknown charge
    %we consider the potential of the unknown charge as it's own position would be [x y z]
    Px = zeros(1, length(M(:, 1)));
    for j = 1:length(M(:, 1))
        Px(1, j) = get_measure([Q; [x y z]], M(j, :)); %compute only on known charge plus the hipotetical unknown charge at the given position
    end
    Px = sum(Px); 
    
    f = sum(abs(P - Px).^2); %error between the real potential and the potential with the unknown charge at [x y z] position
end