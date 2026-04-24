delta = 0.5;
n = 128;
x0 = rand(n) > delta;
imshow(x0);

rng('default')
horseNeigh = [[2 1]; [2 -1]; [-2 -1]; [-2 1]; [1 2];[-1 -2];[-1 2];[1 -2]];
horseNeigh2 = furtherNeighN(5, horseNeigh, horseNeigh);
neigh4 = [[-1 0]; [1 0]; [0 -1]; [0 1]];
neigh8= [[-1 0]; [1 0]; [0 -1]; [0 1]; [1 1]; [-1 -1]; [-1 1]; [1 -1]];
neigh16 = [[0 2]; [0 -2]; [-1 -2]; [-1 2]; [1 -2]; [1 2]; [2 -1]; [-2 -1];
    [-2 1]; [2 1]; [-2 0]; [2 0]; [2 2]; [-2 -2]; [-2 2]; [2 -2]];
neigh8_16 = [neigh16; neigh8];
figure;
hImg = imshow(x0, []);
title('Iteration: 0');
drawnow;
tic
x = BitImgAnneal(x0, 1000, 500000, neigh8, @rotate90Flip, @neighborEnergyDiffColorFcn, hImg);
toc
imshow(x);
function x = BitImgAnneal(x0, initTemp, iter,neighborhood, flipFcn, neighborEnergyFcn, hImg)
    best = x0;
    best_eval = 0;
    current = best;
    current_eval = best_eval;
    n = length(x0);

    for k = 1:iter
        t = initTemp / log(k);
        flippedBits = flipFcn(current);
        flipMask = false(n);
        if ~isempty(flippedBits)
            idx = sub2ind([n n], flippedBits(:,1), flippedBits(:,2));
            flipMask(idx) = true;
        end
        deltaEnergy = neighborFlip(flippedBits, flipMask, current, neighborhood, neighborEnergyFcn);
        if deltaEnergy < 0 || rand() < exp((-deltaEnergy) / t)
            current = flipBits(flippedBits, current);
            
            current_eval = current_eval + deltaEnergy;
            if current_eval < best_eval
                best = current;
                best_eval = current_eval;
            end
        end
        if mod(k, 10000) == 0
            set(hImg, 'CData', current);
            title(['Iteration: ', num2str(k)]);
            drawnow;
            exportgraphics(gcf, "myAnimation.gif", Append=true); % tworzenie gifa
        end
    end
    x = best;
end

function flippedBits = swap2Flip(mat)
    n = length(mat);
    i1 = randi(n, [1 2]);
    i2 = randi(n, [1 2]);
    if mat(i1(1), i1(2)) ~= mat(i2(1), i2(2))
        flippedBits = [i1; i2];
    else
        flippedBits = [];
    end
end

function flippedBits = rotate90Flip(mat)
    n = length(mat);
    l = randi(round(sqrt(n)));
    i1 = randi(n-l);
    i2 = randi(n-l);
    flip_mask = xor(rot90(mat(i1:i1+l, i2:i2+l)), mat(i1:i1+l, i2:i2+l));
    [row, col] = find(flip_mask == 1);
    flippedBits = [row+i1-1 col+i2-1];
end

function flippedBits = transposeFlip(mat)
    n = length(mat);
    l = randi(round(sqrt(n)));
    i1 = randi(n-l);
    i2 = randi(n-l);
    flip_mask = xor(mat(i1:i1+l, i2:i2+l)', mat(i1:i1+l, i2:i2+l));
    [row, col] = find(flip_mask == 1);
    flippedBits = [row+i1-1 col+i2-1];
end

function retmat = flipBits(bits, mat)
    if ~isempty(bits)
        linear_idx = sub2ind(size(mat), bits(:,1), bits(:,2));
        mat(linear_idx) = ~mat(linear_idx);
        retmat = mat;
    else
        retmat = mat;
    end
end

function deltaEnergy = neighborFlip(flips, flipMask, mat, neighborhood_offsets, neighborEnergyFcn)
    n = length(mat);
    deltaEnergy = 0;
    for i = 1:size(flips,1)
        point = flips(i, :);
        for j = 1:size(neighborhood_offsets,1)
            neighbor_offset = neighborhood_offsets(j, :);
            neighbor = point + neighbor_offset; 
            if neighbor(1) < 1 || neighbor(2) < 1 || neighbor(1) > n || neighbor(2) > n
                continue
            end
            deltaEnergy = deltaEnergy + neighborEnergyFcn(point, neighbor, flipMask, mat);
        end
    end
end

function deltaEnergy = neighborEnergySameColorFcn(i1,i2, flipMask, mat)
    deltaEnergy = 0;
    x = i1(1);
    y = i1(2);
    neighbor_x = i2(1);
    neighbor_y = i2(2);
    color = mat(x,y); 
    neighbor_color = mat(neighbor_x, neighbor_y);
    if flipMask(x, y) ~= flipMask(neighbor_x, neighbor_y)
        if color == neighbor_color
            deltaEnergy = deltaEnergy + 1;
        else
            deltaEnergy = deltaEnergy - 1;
        end
    end
end

function deltaEnergy = neighborEnergyDiffColorFcn(i1,i2, flipMask, mat)
    deltaEnergy = 0;
    x = i1(1);
    y = i1(2);
    neighbor_x = i2(1);
    neighbor_y = i2(2);
    color = mat(x,y); 
    neighbor_color = mat(neighbor_x, neighbor_y);
    if flipMask(x, y) ~= flipMask(neighbor_x, neighbor_y)
        if color == neighbor_color
            deltaEnergy = deltaEnergy - 1;
        else
            deltaEnergy = deltaEnergy + 1;
        end
    end
end


function deltaEnergy = neighborEnergyFcn2(i1,i2, flipMask, mat)
    deltaEnergy = 0;
    x = i1(1);
    y = i1(2);
    neighbor_x = i2(1);
    neighbor_y = i2(2);
    color = mat(x,y); 
    neighbor_color = mat(neighbor_x, neighbor_y);
    if flipMask(x, y) ~= flipMask(neighbor_x, neighbor_y)
        if color == neighbor_color
            if norm(i1 - i2) > 5
                deltaEnergy = deltaEnergy + 1;
            else
                deltaEnergy = deltaEnergy - 1;
            end
        else
            if norm(i1 - i2) > 5
                deltaEnergy = deltaEnergy - 1;
            else
                deltaEnergy = deltaEnergy + 1;
            end
        end
    end
end


function deltaEnergy = neighborEnergyFcnZebra(i1,i2, flipMask, mat)
    deltaEnergy = 0;
    x = i1(1);
    y = i1(2);
    neighbor_x = i2(1);
    neighbor_y = i2(2);
    color = mat(x,y); 
    neighbor_color = mat(neighbor_x, neighbor_y);
    if flipMask(x, y) ~= flipMask(neighbor_x, neighbor_y)
        if color == neighbor_color % beda mialy inny kolor
            if i1(2) == i2(2)
                deltaEnergy = +1;
            else
                deltaEnergy = -1;
            end
        else % beda mialy ten sam
            if i1(2) == i2(2)
                deltaEnergy = -1;
            else
                deltaEnergy = +1;
            end
        end
    end
end

function deltaEnergy = neighborEnergyFcnMandel(i1,i2, flipMask, mat)
    n = length(mat);
    deltaEnergy = 0;
    x = i1(1);
    y = i1(2);
    neighbor_x = i2(1);
    neighbor_y = i2(2);
    color = mat(x,y); 
    neighbor_color = mat(neighbor_x, neighbor_y);
    [mand_x, mand_y] = map_to_mandelspace(x, y, n);
    [mand_neighbor_x, mand_neighbor_y] = map_to_mandelspace(neighbor_x, neighbor_y, n);

    z = 0;
    i = 0;
    neighbor_diverges =false;
    p_diverges = false;
    while i < 100
        z = z*z + (mand_x + mand_y * 1i);
        if abs(z) > 2
            p_diverges = true;
            break
        end
        i = i + 1;
    end
    z = 0;
    i = 0;
    while i < 100
        z = z*z + (mand_neighbor_x + mand_neighbor_y * 1i);
        if abs(z) > 2
            neighbor_diverges = true;
            break
        end
        i = i + 1;
    end

    if flipMask(x, y) ~= flipMask(neighbor_x, neighbor_y) 
        if color == neighbor_color % beda mialy inny kolor
            if xor(neighbor_diverges, p_diverges) % sa po innych stronach
                if (p_diverges && color)
                    deltaEnergy = -1;
                else
                    deltaEnergy = +1;
                end
            else % sa po tej samej stronie
                deltaEnergy = 0;
            end
        else % beda mialy ten sam
            if xor(neighbor_diverges, p_diverges) % sa po innych stronach
                deltaEnergy = 0;
            else % sa po tej samej stronie
                if p_diverges && color
                    deltaEnergy = -1;
                else
                    deltaEnergy = +1;
                end
            end
        end
    else

    end

    function [x, y] = map_to_mandelspace(x, y, n)
        x = (2.5 / n) * x - 2;
        y = (2.4 / n) * y - 1.2;
    end
end

function out = furtherNeighN(N, recursive_neigh, og_neigh)
    if N == 0
        out = unique(recursive_neigh, "rows");
        return
    end
    out = recursive_neigh;
    for i = 1:length(og_neigh)
        out = [out; recursive_neigh + og_neigh(i, :)];
    end
    out = furtherNeighN(N - 1, out, og_neigh);
end
