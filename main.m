import Simplex

close all

%initializing Simplex class
s = Simplex(@cost, {@bound1}, 10, [-17 -17 10], 1e-10, 150);
s.dt = 0.1; %animation delta time between frames (0 = off)
s.field = 20; %figure subspace of view
s.slices = 20; %15 planes to draw the isolevel maps
s.color = 'magenta'; %simplex politope color
s.plot = true; %enable the polotting

[value, coordinates, flips, halvings, area] =  s.compute() %compute the algorithm

function f = bound1(x, y, z) %sphere bound
    f = -((x+15).^2 + (y+15).^2 + (z-10).^2 - 20^2);
end

function f = get_measure(Q, m) %get the potential from the measure in a certain position
    m = repmat(m, [length(Q(1:end, :)), 1]); %adapt the matrix
    D = sqrt(sum((Q-m).^2, 2)); %norm
    
    f = 1/sum(D); %potential
end

function f = cost(x, y, z)
    Q = [1 0 0; 0 1 0; 0 0 1; 0 0 0; 0 0 1; 0 1 1; 1 0 0; 1 1 1;]; %known charges
    q = [0.5 0.5 0.5]; %unknown charge
    M = 2*[2 0 0; 0 2 0; 0 0 2]; %measures coordinates in the R^3 space
    
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