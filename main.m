import Simplex

close all
clear
clc

global z
global M
global Q
global q
global C
global c

z  = -.2;
C = [.1 .3 -.4 .3 -4. .1 .2 -.2 .1 .4];
Q = [-.4 -.4 -.4
    0 -.2 -.3
    .4 0 .3
    .3 .3 .3
    -.2 .3 -.2
    .2 -.2 -.1
    -.2 .1 0
    -.3 .4 .4
    .3 .4 .3
    -.3 .4 z]; %10 charges
q = Q(end, :); %the unknown charge
Q = Q(1:end-1, :); %the other 9 known charges
c = C(end);
C = C(1:end-1);
M = [];

%mapping the measurers box
for i = -1:1:1
    for j = -1:1:1
        %X planes
        M = [M; [-1, i, j]];
        M = [M; [1, i, j]];

        %Y planes
        M = [M; [i, -1, j]];
        M = [M; [i, 1, j]];

        %Z planes
        M = [M; [i, j, -1]];
        M = [M; [i, j, 1]];
    end
end

%initializing Simplex class
s = Simplex(@cost_function, {}, [-.2 -.2 .1], .075, 1e-5, 25); % NOTE: I have changed the value
s.dt = 0; %animation delta time between frames (0 = off)
s.field = 1; %figure subspace of view
s.slices = 10; %15 planes to draw the isolevel maps
s.color = 'green'; %simplex politope color
s.plot = true; %enable the polotting

[value, minimum, flips, halvings, area] =  s.compute() %compute the algorithm


disp('known charge');   
disp(q);

error_v = abs(cost_function(minimum) - cost_function([q(1) q(2), c]));
error_c = abs(minimum - [q(1) q(2), c]);
disp('error');
disp(error_v)
disp(error_c)



function f = bound1(x, y, z) %sphere bound
    f = -((x+.6).^2 + (y-.6).^2 + (z-.1).^2 - 1^2); %sphere
end

function f = bound2(x, y, z) %sphere bound
    f = -((x+1.5).^2 + (y+1.5).^2 + (z+.5).^2 - 2^2);
end

function E = cost_function(qu)
    global M 
    global Q 
    global q
    global C
    global c
    global z
    
    U = get_all_potentials([Q; q], [C c]); %get the whole potentials measured from the measurers
    Un = get_all_potentials(Q, C); %get the known charges potentials
    Uc = U - Un; %Uc get only the potential of the unknown charge
    %In this way we'll have only the difference from the unknown charge's
    %potential and the generic charge computed in the qu position
    %We don't use the absolute value because in a further step we'll
    %exponentiate by 2 all the values
    
    Uc = abs(Uc - get_all_potentials([qu(1), qu(2), z], qu(3)));
    Uc = Uc./(length(M)*max(C));
    
    E = sum(Uc.^2); %sum quadratically the single cost functions to minimize them all
end

function U = get_all_potentials(Q, C)
    global M
    
    U = [];
    for i = 1:length(M(:, 1)) %for each measurer
        Um = 0;
        for j = 1:length(Q(:, 1)) %for each charge
            Um = Um + get_potential(Q(j, :), C(j), M(i, :)); %sum the potential of the j-th charge
        end
        
        U = [U Um]; %store the measurer detection
    end
end

function u = get_potential(p, c, m) %compute charge's potential detected from a measurer
    u = c/(norm(p-m)); %yet normalized
end