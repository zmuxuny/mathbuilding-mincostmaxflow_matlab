
```matlab
clear;
clc;
V = [
     0 5 8 0 0
     0 0 0 3 4
     0 2 0 10 0
     0 0 0 0 8
     0 0 0 0 0
     ]; %volume
C = [
     0 8 7 0 0
     0 0 0 2 9
     0 5 0 9 0
     0 0 0 0 4
     0 0 0 0 0
     ]; %cost
n = size(V, 2);

flow = 0;
cost = 0;

f = zeros(n, n); % 初始流量

% 创建图
G = digraph(V);
highlightedEdges = []; % 用于保存最小费用流的路径
a = inf * (ones(n, n) - eye(n, n)); % 有向加权图
% 这行代码创建了一个有向加权图的邻接矩阵。在这个矩阵中，每个元素代表了对应的两个节点之间的边的权重。
%onepre(n, n)创建了一个n行n列的矩阵，所有元素都是1。eye(n, n)创建了一个n行n列的单位矩阵，对角线上的元素是1，其他元素都是0。
%onepre(n, n) - eye(n, n)的结果是一个n行n列的矩阵，对角线上的元素是0，其他元素都是1。
%然后，这个矩阵被乘以inf，所以最终的结果是一个n行n列的矩阵，对角线上的元素是0，
%其他元素都是inf。这表示在初始状态下，所有的节点都是不连通的，除了每个节点自身（对角线上的元素）。

% !用费用作为边权
for i = 1:n

    for j = 1:n

        a(i, j) = C(i, j);
        a(j, i) = -C(i, j);
    end

end

while 1

    %spfa找最短路
    p = inf * ones(1, n); % 初始访问向量 %!表示的是从源点到某个点的路径长度，初始化为inf
    p(1) = 0; %1号点是源点，它到自己距离为0
    pre = 1:n; %s记录前驱节点
    pre(1) = 0;
    mf = zeros(1, n); %记录容量上限
    mf(1) = inf;
    q = zeros(1, 100000); %队列
    hh = 1;
    tt = 1;
    vis = zeros(1, 10000); %记录有没有入队
    q(tt) = 1;
    tt = tt + 1;
    vis(1) = 1;

    while (hh <= tt)
        u = q(hh);
        hh = hh + 1;

        if isnumeric(u) && u > 0 && round(u) == u
            vis(u) = 0;
        else
            break;
        end

        for j = 1:n

            if j == u || j == 1
                continue;
            end

            if V(u, j) > 0

                if p(j) > p(u) + a(u, j)
                    p(j) = p(u) + a(u, j);
                    pre(j) = u;

                    if mf(u) <= V(u, j)
                        mf(j) = mf(u);
                    else
                        mf(j) = V(u, j);
                    end

                    if vis(j) ~= 1
                        q(tt) = j;
                        tt = tt + 1;
                        vis(j) = 1;
                    end

                end

            end

        end

    end

    if mf(n) == 0
        disp("找不到增广路了");
        break;
    end %表示没有从源点1到汇点n的路径

    v = n;
    disp("当前的mf(n)");
    mf(n)

    while v ~= 1
        u = pre(v);
        V(u, v) = V(u, v) - mf(n);
        f(u, v) = f(u, v) + mf(n);
        V(v, u) = V(v, u) + mf(n);
        v = u;
        flow = flow + mf(n);
        cost = mf(n) * p(n) + cost;
    end

    figure

    path = [];
    t = n;

    while t ~= 1
        path = [pre(t), t; path];
        t = pre(t);
    end

    highlightedEdges = [highlightedEdges; path];

    % 创建流量和费用的字符串数组
    edgeLabels = strings(size(G.Edges, 1), 1);

    for i = 1:size(G.Edges, 1)
        edgeLabels(i) = sprintf('%d, %d', f(G.Edges.EndNodes(i, 1), G.Edges.EndNodes(i, 2)), C(G.Edges.EndNodes(i, 1), G.Edges.EndNodes(i, 2)));
    end

    % 使用新的字符串数组作为'EdgeLabel'的值
    h = plot(G, 'EdgeLabel', edgeLabels);

    highlight(h, highlightedEdges(:, 1), highlightedEdges(:, 2), 'EdgeColor', 'r'); % 将最小费用流的路径标红
    title('流量分布')

    drawnow;
end

f % 流量矩阵
flow % 最大流量值
cost % 最小流量值

```
