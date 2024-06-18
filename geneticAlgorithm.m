% 遺伝的アルゴリズムのシミュレーション

function geneticAlgorithm
    populationSize = 40;    % 個体数
    generations = 100;  % 世代数 (ループ数)
    crossoverRate = 0.7;    % 交叉率
    mutationRate = 0.05;    % 突然変異率
    tournamentSize = 3; % トーナメントのサイズ
    gifFilename = 'genetic_algorithm_three_hump_animation.gif';
    optimalValue = 0; % 目的関数の最適値
    % seed = 1;   % 局所解
    % seed = 2;   % global
    seed = 810;

    rng(seed);  % seed値を固定
    
    % 初期集団(n=40)をランダムに[-5, 5]の範囲で生成
    population = 10 * rand(populationSize, 2) - 5; % [-5, 5]の範囲

    % 等高線図の準備
    [X, Y] = meshgrid(-2:0.01:2, -2:0.01:2);
    Z = myfunc(X, Y);

    % 各世代の最良適応度を記録する配列
    bestFitnessHistory = zeros(generations, 1);
    errorHistory = zeros(generations, 1);

    figure;
    for gen = 1:generations

        % 個体数のチェック
        % disp(['Generation: ' num2str(gen) ', Population Size: ' num2str(size(population, 1))]);

        % 適応度の計算
        fitness = evaluateFitness(population);

        % 個体群の進化
        % selectedPopulation = selection(population, fitness);
        selectedPopulation = selection(population, fitness, tournamentSize);
        crossedPopulation = crossover(selectedPopulation, crossoverRate);
        mutatedPopulation = mutation(crossedPopulation, mutationRate);
        population = mutatedPopulation;

        % 最良適応度の更新
        bestFitnessHistory(gen) = max(fitness);
        errorHistory(gen) = abs(optimalValue - bestFitnessHistory(gen));

        % プロット
        subplot(2, 1, 1); % 上部の等高線図
        contour(X, Y, Z, 20);
        hold on;
        scatter(population(:, 1), population(:, 2), 'ro');

        title(['Generation: ' num2str(gen)]);
        xlim([-2, 2]);
        ylim([-2, 2]);
        hold off;

        subplot(2, 1, 2); % 下部の誤差の散布図
        scatter(1:gen, errorHistory(1:gen), 'ro');
        xlim([0, generations]);
        ylim([0 max(errorHistory)]);
        xlabel('Generation');
        ylabel('Error from Optimal Value');
        title('Error of Best Fitness Over Generations');
        drawnow;

        % GIFのフレームを保存
        frame = getframe(gcf);
        im = frame2im(frame);
        [imind, cm] = rgb2ind(im, 256);

        if gen == 1
            imwrite(imind, cm, gifFilename, 'gif', 'Loopcount', inf, 'DelayTime', 0.2);
        else
            imwrite(imind, cm, gifFilename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.2);
        end
    end

    % 最終的な適応度履歴のプロット
    figure;
    plot(bestFitnessHistory, 'r-', 'LineWidth', 2);
    xlabel('Generation');
    ylabel('Best Fitness');
    title('Best Fitness over Generations');
end

function z = myfunc(x, y)
    z = 2 * x.^2 - 1.05 * x.^4 + x.^6 / 6 + x .* y + y.^2;    % Three-Hump Camel Function
    
    % z = x.^2 + y.^2;

    % ackley (0 at x=0)
    % z = -20 * exp(-0.2 * sqrt(0.5 * (x.^2 + y.^2))) - exp(0.5 * (cos(2 * pi * x) + cos(2 * pi * y))) + exp(1) + 20;
    
    %rastrigin (0 at x=0)
    % A = 10;
    % z = A * 2 + (x.^2 - A * cos(2 * pi * x)) + (y.^2 - A * cos(2 * pi * y));
end

% 適応度の計算
function fitness = evaluateFitness(population)
    fitness = -myfunc(population(:, 1), population(:, 2));  % 目的関数の負
end

% 選択処理
% function selectedPopulation = selection(population, fitness)
%     populationSize = size(population, 1);
%     selectedPopulation = zeros(size(population));
%     totalFitness = sum(fitness);
% 
%     for i = 1:populationSize
%         % ルーレット選択
%         pick = rand * totalFitness;
%         current = 0;
%         for j = 1:populationSize
%             current = current + fitness(j);
%             if current > pick
%                 selectedPopulation(i, :) = population(j, :);
%                 break;
%             end
%         end
%     end
% end

% 個体の選択処理 (トーナメント選択)
function selectedPopulation = selection(population, fitness, tournamentSize)
    populationSize = size(population, 1); % 集団のサイズ（個体数）
    selectedPopulation = zeros(size(population)); % 選択された個体を格納する配列の初期化

    for i = 1:populationSize
        % トーナメントに参加する個体をランダムに選択
        competitorsIndex = randi(populationSize, [tournamentSize, 1]); % トーナメント参加者のインデックスをランダムに選ぶ
        competitorsFitness = fitness(competitorsIndex); % トーナメント参加者の適応度を取得

        % 最も適応度の高い個体を選択
        [~, bestIndex] = max(competitorsFitness); % トーナメント参加者の中で最も適応度が高い個体のインデックスを取得
        selectedPopulation(i, :) = population(competitorsIndex(bestIndex), :); % 最も適応度の高い個体を選択された集団に追加
    end
end

% 交叉処理
function crossedPopulation = crossover(selectedPopulation, crossoverRate)
    populationSize = size(selectedPopulation, 1); % 選択された集団のサイズ（個体数）
    crossedPopulation = selectedPopulation; % 交叉後の集団を選択された集団で初期化

    for i = 1:2:populationSize
        if rand < crossoverRate % 交叉率に基づいて交叉を行うか決定
            crossoverPoint = randi([1 size(selectedPopulation, 2)-1]); % 交叉点をランダムに決定
            % 交叉点以降の部分を入れ替える
            temp = crossedPopulation(i, crossoverPoint+1:end); % 親1の交叉点以降の部分を一時的に保存
            crossedPopulation(i, crossoverPoint+1:end) = crossedPopulation(i+1, crossoverPoint+1:end); % 親2の交叉点以降の部分を親1にコピー
            crossedPopulation(i+1, crossoverPoint+1:end) = temp; % 親1の交叉点以降の部分を親2にコピー
        end
    end
end

% 突然変異処理
function mutatedPopulation = mutation(crossedPopulation, mutationRate)
    populationSize = size(crossedPopulation, 1); % 交叉後の集団のサイズ（個体数）
    mutatedPopulation = crossedPopulation; % 突然変異後の集団を交叉後の集団で初期化

    for i = 1:populationSize
        for j = 1:size(crossedPopulation, 2)
            if rand < mutationRate % 突然変異率に基づいて突然変異を行うか決定
                % 突然変異により現在の値に小さなランダムな値を加える
                mutatedPopulation(i,j) = mutatedPopulation(i,j) + (rand - 0.5) * 0.2; % -0.1から0.1の範囲で変化を加える
            end
        end
    end
end

