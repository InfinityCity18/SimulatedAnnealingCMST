rng('default')
total = 0;
k = 1000;
for i = 1:k
    sudoku1 = getSudokuAtLine("sudoku_long.txt", randi(800000));
    fixedMask = (sudoku1 ~= 0);
    sudoku1 = SudokuX0Squares(sudoku1, fixedMask);
    % options = optimoptions( ...
    %     "simulannealbnd", ...
    %     "PlotInterval", 1000, ... % Wybierz interwał wizualizacji
    %     "ReannealInterval", 0, ...
    %     "AnnealingFcn", @(optimValues,~) SudokuNeighborSquareSwap2(optimValues.x, fixedMask), ... % Wybierz funkcję generującą sąsiada
    %     "TemperatureFcn","temperaturefast", ...
    %     "ObjectiveLimit", 1, ...
    %     "DataType","custom", ...
    %     'InitialTemperature', 20000, ... % Wybierz temperaturę startową
    %     "MaxStallIterations", 100000, ... % Maks iteracji
    %     "FunctionTolerance", 1e-15, ...% Tolerancja (wyżarzanie zakończy się gdy średnia zmiana wartości funkcji będzie mniejsza od tolerancji)
    %     "PlotFcn", {@(~, optimvalues, state) plotSudoku(0, optimvalues.x, state, fixedMask), @saplotf, @saplotbestf} ...
    % );
    tic
    %[x, fval] = simulannealbnd(@FastSudokuObjective, sudoku1, [], [], options);
    x = SudokuAnnealOwn(@FastSudokuObjective, @SudokuNeighborSquareSwap2, fixedMask, sudoku1, 2200, 300000, 7000);
    disp(FastSudokuObjective(x));
    fprintf("i = %d\n", i);
    if FastSudokuObjective(x) == 0
        total = total + 1;
    end
    toc
end
disp(total);

function fval = SudokuObjective(sudoku)
    fval = 0;
    for row = 1:3:9
        for col = 1:3:9
            V = sudoku(row:row+2, col:col+2);
            [~,~,ix] = unique(V);
            fval = fval + sum(accumarray(ix,1) - 1);
        end
    end
    for i = 1:9 % kolumny i wiersze
        V = sudoku(1:9, i);
        [~,~,ix] = unique(V);
        fval = fval + sum(accumarray(ix,1) - 1);
        V = sudoku(i, 1:9);
        [~,~,ix] = unique(V);
        fval = fval + sum(accumarray(ix,1) - 1);
    end
end

function x = SudokuAnnealOwn(objectiveFcn, neighborFcn, fixedMask, x0, initTemp, max_iter, renneal_threshold)
    best = x0;
    best_eval = objectiveFcn(x0);
    current = best;
    current_eval = best_eval;
    plot_x = [];
    plot_y = [];
    reanneal_counter = 0;
    k = 0;
    iter = 0;

    while iter < max_iter
        t = initTemp / (k + 1);
        candidate = neighborFcn(current,fixedMask);
        candidate_eval = objectiveFcn(candidate);
        if candidate_eval < best_eval || rand() < exp((current_eval - candidate_eval) / t)
            current = candidate;
            current_eval = candidate_eval;
            if candidate_eval < best_eval
                reanneal_counter = 0;
                best = candidate;
                best_eval = candidate_eval;
                if best_eval == 0
                    x = best;
                    return
                end
            else 
                reanneal_counter = reanneal_counter + 1;
            end
        end
        if reanneal_counter > renneal_threshold
            k = 0;
            reanneal_counter = 0;
        end
        if false && mod(iter, 10000) == 0
            tiledlayout(2,1);
            plot_x = [plot_x iter];
            plot_y = [plot_y current_eval];
            nexttile
            plot(plot_x,plot_y);
            nexttile
            plotSudoku(0, current, 'iter', fixedMask);
        end
        k = k + 1;
        iter = iter + 1;
    end
    x = best;
end

function fval = FastSudokuObjective(sudoku)
    % Since our Neighbor function ensures each 3x3 block is a permutation of 1-9,
    % we only need to check for duplicates in rows and columns.
    
    % Create a 9x9x9 logical cube where (i,j,k) is true if sudoku(i,j) == k
    % This is extremely fast in MATLAB
    vals = reshape(1:9, 1, 1, 9);
    matches = (sudoku == vals);
    
    % Row penalty: Count how many times each number is missing from each row
    % If a row has all 1-9, sum(matches, 2) will be all 1s.
    row_counts = sum(matches, 2); 
    row_penalty = sum(row_counts == 0, 'all');
    
    % Column penalty
    col_counts = sum(matches, 1);
    col_penalty = sum(col_counts == 0, 'all');
    
    fval = row_penalty + col_penalty;
end

function sudoku = SudokuNeighborRowSwap2(sudoku, fixedMask)
    row_i = randi(9);
    i1 = randi(9);
    i2 = randi(9);
    while fixedMask(row_i, i1) == 1 || fixedMask(row_i, i2) == 1 % oba nie sa fixed
        i1 = randi(9);
        i2 = randi(9);
    end
    temp = sudoku(row_i, i1);
    sudoku(row_i, i1) = sudoku(row_i, i2);
    sudoku(row_i, i2) = temp;
end

function sudoku = SudokuNeighborSquareSwap2(sudoku, fixedMask)
    choices = [1 4 7];
    row_i = choices(randi(3));
    col_i = choices(randi(3));
    i1 = randi(9);
    i2 = randi(9);
    subsquare = sudoku(row_i:row_i+2, col_i:col_i+2);
    subFixedMask = fixedMask(row_i:row_i+2, col_i:col_i+2);
    while subFixedMask(i1) == 1 || subFixedMask(i2) == 1 || sum(subFixedMask, 'all') >= 8  % jeden jest fixed lub >= 8 jest fixed
        row_i = choices(randi(3));
        col_i = choices(randi(3));
        subsquare = sudoku(row_i:row_i+2, col_i:col_i+2);
        subFixedMask = fixedMask(row_i:row_i+2, col_i:col_i+2);
        i1 = randi(9);
        i2 = randi(9);
    end
    subsquare([i1 i2]) = subsquare([i2 i1]);
    sudoku(row_i:row_i+2, col_i:col_i+2) = subsquare;
end

function sudoku = SudokuX0Row(sudoku, fixedMask)
    for i = 1:9
        possible_choice = setdiff(1:9, sudoku(i, :));
        after_shuffling = possible_choice(randperm(length(possible_choice)));
        row = sudoku(i, :);
        row(~fixedMask(i, :)) = after_shuffling;
        sudoku(i, :) = row;
    end
end

function sudoku = SudokuX0Squares(sudoku, fixedMask)
    for i = 1:3:9
        for j = 1:3:9
            fixedMaskSquare = fixedMask(i:i+2, j:j+2);
            possible_choice = setdiff(1:9, sudoku(i:i+2, j:j+2));
            after_shuffling = possible_choice(randperm(length(possible_choice)));
            square = sudoku(i:i+2, j:j+2);
            square(~fixedMaskSquare) = after_shuffling;
            sudoku(i:i+2, j:j+2) = square;
        end
    end
end
        

function [grid, difficulty] = getSudokuAtLine(filePath, n)
    % filePath: Path to the .txt file
    % n: The specific line number to retrieve
    
    fid = fopen(filePath, 'r');
    if fid == -1, error('File not found.'); end
    
    % Ensure the file closes even if the function is interrupted
    cleanup = onCleanup(@() fclose(fid));

    % Skip the first n-1 lines
    for i = 1:(n-1)
        if ~ischar(fgetl(fid))
            error('The file ended before reaching line %d.', n);
        end
    end

    % Read the target line
    line = fgetl(fid);
    if ~ischar(line), error('Line %d is empty or doesn''t exist.', n); end
    
    % Split the string and difficulty
    parts = strsplit(strtrim(line), ' ');
    sudokuStr = parts{1};
    
    % Set difficulty (default to NaN if missing)
    difficulty = NaN;
    if numel(parts) > 1
        difficulty = str2double(parts{2});
    end

    % Replace 'x' with '0' and convert to numeric 9x9
    % We do this character-wise to stay fast
    sudokuStr(sudokuStr == 'x') = '0';
    grid = reshape(sudokuStr - '0', 9, 9)';
end

function stop = plotSudoku(~, x, state, fixedMask)
    stop = false; % Set to true if you want to terminate the solver
    
    % Current board state (x is usually a vector or matrix depending on setup)
    currentX = x;
    if isvector(currentX)
        board = reshape(currentX, 9, 9);
    else
        board = currentX;
    end

    switch state
        case 'iter'
            hold on;
            cla;
            
            % Setup axis
            axis equal;
            set(gca, 'YDir', 'reverse', 'XLim', [0.5 9.5], 'YLim', [0.5 9.5]);
            
            % 1. Draw Grid Lines
            for i = 0.5:1:9.5
                lw = 1;
                if mod(i-0.5, 3) == 0, lw = 3; end
                line([0.5 9.5], [i i], 'Color', 'k', 'LineWidth', lw);
                line([i i], [0.5 9.5], 'Color', 'k', 'LineWidth', lw);
            end

            % 2. Plot Numbers
            for r = 1:9
                for c = 1:9
                    val = board(r, c);
                    if val ~= 0
                        % Visual distinction for initial numbers
                        if fixedMask(r, c)
                            col = 'r'; fw = 'bold';
                        else
                            col = 'white'; fw = 'normal';
                        end
                        text(c, r, num2str(val), 'HorizontalAlignment', 'center', ...
                            'FontSize', 12, 'FontWeight', fw, 'Color', col);
                    end
                end
            end
            drawnow limitrate;
            
        case 'done'
            hold off;
    end
end