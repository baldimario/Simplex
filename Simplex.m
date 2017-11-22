classdef Simplex
    properties
        func = 0; %function to minimize
        bounds = {}; %bounds list
        start_area = 5; %triangle start area
        start_point = [0 0]; %triangle start point
        stopping_area = 1e-5; %first stop condition on the minimum area
        max_steps = 20; %second stop condition, maximum halvings count
        dt = 0; %animation delta time between frames (0 = off)
        plot = false; %enable or disable the plotting
        field = 10; %figure subspace of view
        slices = 5; %number of planes to draw the isolevel maps
        color = 'blue'; %simplex polytope color
    end
    methods
        %contructor
        function obj = Simplex(func, bounds, start_area, start_point, stopping_area, max_steps)
            obj.func = func;
            obj.bounds = bounds;
            obj.start_area = start_area;
            obj.start_point = start_point;
            obj.stopping_area = stopping_area;
            obj.max_steps = max_steps;
        end
        
        %core function
        function [value, coordinates, flips, halvings, area] = compute(obj)
            Polytopes = []; %polytopes list
            halvings = 1; %halving counter
            flips = 1; %flip counter
            persist = true; %loop condition

            P = obj.get_first_polytope(); %getting the fist polytope
            [Polytopes, flips] = obj.add_polytope(Polytopes, flips, P);
            
            if obj.plot %if plot it's enables 
                %let's draw
                figure; 
                hold on;
                
                obj.draw_function();
                obj.draw_polytope(P);
                obj.draw_bounds();
                view(65,8);
            end
            
            while persist
                m = obj.find_maximums(P); %find the maximums positions to select the right vertex
                for j = 1:4 %for each vertex
                    vx = obj.predict_vertex(P, m(j)); %predict the flipped coordinates

                    if obj.check_bounds(vx) %check if the prediction respect the bounds
                        P = obj.flip_polytope(P, m(j)); %flip the right vertex
                        break;
                    end
                end

                [Polytopes, flips] = obj.add_polytope(Polytopes, flips, P); %add the new polytope to the list
                
                if obj.plot %if plot enabled draw the polytope
                    obj.draw_polytope(P);
                end
                
                if(obj.dt ~= 0 && obj.plot) %if animation it's enabled (dt != 0) (if plot is disabled then disable animations)
                    pause(obj.dt); %pause for dt seconds
                end

                if obj.watchdog(Polytopes) %check if the last polytopes are are in a flipping loop condition
                    P = obj.halve(P, j);  %halve the polytope on the last vertex pivot
                    halvings = halvings + 1; %increment the halvings counter
                end
                
                if flips >= obj.max_steps %if the maximum halvings counter reach the limit
                    persist = false; %break the loop
                    continue;
                end
                
                persist = ~obj.detect_stop(P); %if the area reach the stop area condition break the loop
            end

            area = obj.get_area(P); %return the area
            [value, coordinates] = obj.get_results(P); %return the value of the minimum and the coordinates
        end

        %predict_vertex gets the polytope and the vertex index to flip and
        %returns the coordinates of the flipped vertex
        function vx = predict_vertex(obj, P, j)
            vx = P(j, 1:end); %get the pivot
            
            %get other vertex
            V = P([1:j-1 j+1:end], 1:end); 
            v1 = V(1, 1:end);
            v2 = V(2, 1:end);
            v3 = V(3, 1:end);
            
            
            p = (v1 + v2 + v3)/3; %find the middle point
            d = (p-vx); %get the distance from the pivot and the middle point
            vx = p+d; %flip the vertex in direction of middle point
        end

        %draw the bounds in according to the field and slices parameters
        function y = draw_bounds(obj)
            if(length(obj.bounds) == 0)
                return
            end
            
            X = -obj.field:obj.field/obj.slices:obj.field;
            Y = -obj.field:obj.field/obj.slices:obj.field;
            Z = -obj.field:obj.field/obj.slices:obj.field;
            V = zeros(length(obj.bounds), obj.slices, length(X), length(Y));
            
            for k = 1:length(Z)
                J = zeros(length(X), length(Y));
                
                for b = 1:length(obj.bounds)
                    B = zeros(length(X), length(Y));
                    
                    for i = 1:length(X)
                        for j = 1:length(Y)
                            B(i, j) = obj.bounds{b}(X(i), Y(j), Z(k));
                        end
                    end
                    Bm = B >= 0;
                    Bm = double(Bm);
                    B(Bm == 0) = NaN;
                    J = J.*B;
                end
                
                offset = ((k-(length(Z)/2))*(obj.field/obj.slices));
                K = (ones(length(X), length(Y))*offset)+(obj.field/obj.slices)/2;
                
                
                colormap(jet);
                hold on
                s = surf(X, Y, K, 'FaceAlpha', 0.5,'LineStyle','none', 'cdata', J, 'cdatamapping', 'direct'); 
                colormap(jet);
                hold on
                colormap(jet);
            end
        end

        %get_first_polytope, as the name shows, computes the first polytope
        %according to the given initial point
        function P = get_first_polytope(obj)
            l = sqrt(4*obj.start_area/sqrt(3)); %get the side of the equilateral triangle
            h = sqrt(3)*l/2; %get the height of the equilateral triangle

            v1 = [obj.start_point]; %the first vertex is the start point
            v2 = [obj.start_point + [h, 0 h]]; 
            v3 = [obj.start_point + [h, l/2 0]];
            v4 = [obj.start_point + [h, -l/2 0]];

            P = [v1; v2; v3; v4];
        end

        %get_results returns the value of the minimum of func calculated on
        %each vertex of the last polytope and its coordinate
        function [z, x] = get_results(obj, P)
            v = obj.func(P(1:end, 1), P(1:end, 2), P(1:end, 3)); %calculate the function for each vertex
            [z, j] = min(v); %find the minimum
            x = P(j, 1:end);
        end

        %halve halves the polytope taking the index of the pivot vertex
        function P = halve(obj, P, j)
            vx = P(j, 1:end); %get the pivot vertex
            
            %get the other vertex
            V = P([1:j-1 j+1:end], 1:end); 
            v1 = V(1, 1:end);
            v2 = V(2, 1:end);
            v3 = V(3, 1:end);

            %compute the middle points of each vertex with the pivot vertex
            v1 = (vx+v1)/2;
            v2 = (vx+v2)/2;
            v3 = (vx+v3)/2;

            P = [vx; v1; v2; v3];
        end

        %get_area returns the given polytope area
        function y = get_area(obj, P)
            l = norm(P(1, 1:end) - P(2, 1:end)); %calculates the side
            y = 4*sqrt(3)*l/2; %single equilateral triangle area multiplied by 4 faces
        end

        %detect_stop returns a boolean value
        %true if the given polytope area it's minus of the stopping area
        function y = detect_stop(obj, P)
            a = obj.get_area(P); %getting the polytope area

            y = a <= obj.stopping_area; %stopping condition
        end

        %add_polytope adds the polytope to the list and increments the flips counter
        function [Polytopes, flips] = add_polytope(obj, Polytopes, flips, P)
            Polytopes(:, :, flips) = P; %add the polytope to the list
            flips = flips+1; %increment the halvings
        end

        %whatchdog gets the polytopes list and check if there is a repeated
        %polytope in the last 8 elements of the list, if so returns true
        %otherwise returns false
        function y = watchdog(obj, Polytopes)
            y =false;
            
            if(length(Polytopes) >= 8)
                for j = 1:8
                    for k = 1:8
                        if j ~= k
                            if Polytopes(:,:,end-8+j) == Polytopes(:,:,end-8+k)
                                y = true;
                            end
                        end
                    end
                end
            end
        end

        %draw_function draws the function to minimize in according to the field and slices parameters
        function y = draw_function(obj)
            X = -obj.field:obj.field/obj.slices:obj.field;
            Y = -obj.field:obj.field/obj.slices:obj.field;
            Z = -obj.field:obj.field/obj.slices:obj.field;
            V = zeros(obj.slices, length(X), length(Y));
            
            %calculates the function values in a 3-dimensional grid
            for k = 1:length(Z)
                for i = 1:length(X)
                    for j = 1:length(Y)
                        V(k, i, j) = -obj.func(X(i), Y(j), Z(k));
                    end
                end
            end
            
            
            for j = 1:length(Z)
                %calculate the plane
                offset = ((j-(length(Z)/2))*(obj.field/obj.slices));
                Z = ones(length(X), length(Y))*offset;
                
                %draw the plane with the iso-level colors of the function
                %computed in its coordinates
                colormap(jet);
                surf(X, Y, Z, reshape(V(j,:,:), [length(X), length(Y)]), 'FaceAlpha', 0.25,'LineStyle','none', 'FaceColor','interp');
            end
        end

        %find_maximums
        function m = find_maximums(obj, P)
            m = [];
            
            %repeat for a number of times equals to polytope vertex
            for j = 1:length(P(:, 1))
                %create a list containing the every polytope vertex
                V = [obj.func(P(1, 1), P(1, 2), P(1, 3)),
                    obj.func(P(2, 1), P(2, 2), P(2, 3)),
                    obj.func(P(3, 1), P(3, 2), P(3, 3)),
                    obj.func(P(4, 1), P(4, 2), P(4, 3))];
                
                %find the maximum
                [r, y] = max(V);
                
                %remove the vertex from the polytope
                P(y, 1:end) = NaN;
                
                %save the maximum
                m = [m y];
            end
        end
        
        %flip_polytope flips the j-th vertex of the polytope 
        function P = flip_polytope(obj, P, j)
            vx = P(j, 1:end); %get the vertex to flip
            
            %get other vertex
            V = P([1:j-1 j+1:end], 1:end); 
            v1 = V(1, 1:end);
            v2 = V(2, 1:end);
            v3 = V(3, 1:end);

            p = (v1+v2+v3)/3; %compute the midde point of other vertex
            d = (p-vx); %calculate the distance vector from the middle point and the vertex to flip
            vx = p+d; %flip the vertex

            P(j, 1:end) = vx;
        end

        %check_bounds returns true if every bound it's satisfied by the v coordinates as their inputs
        function y = check_bounds(obj, v)
            y = true;

            for j = 1:length(obj.bounds)
                if obj.bounds{j}(v(1), v(2), v(3)) < 0
                    y = false;
                end
            end
        end

        %draw_polytope draws the polytope in the 3-d space
        function y = draw_polytope(obj, P)
            X = P(1:end, 1);
            Y = P(1:end, 2);
            Z = P(1:end, 3);

            f = [1 2 3; 1 2 4; 1 3 4; 2 3 4]; %faces list every index it's a vertex of the polytope and 3 vertex form a triangular face of the polytope

            patch('Faces',f,'Vertices',P, 'EdgeColor','black','FaceColor',obj.color,'LineWidth',2, 'FaceAlpha', 0.7);

            y = 1;
        end
    end
end
