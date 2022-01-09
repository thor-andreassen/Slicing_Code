function [free_edges,watertight] = findFreeEdges(edges)
        % Function to determine free edges of mesh and if mesh is
        % watertight
        % function takes in a matrix of edge connections as a (n x 2)
        % matrix, and returns a list of (m x 2) free edges if there are
        % any. Free edges are defined as any edges that are only connected
        % to a single face, and therefore only appear once in the edges
        % list. The function also returns a "watertight" measure, that
        % defines if the entire list of edges is watertight. This is true
        % if there are no free edges, and false if there are free edges.
        sort_edges=sort(edges,2);
        [B,ix] = sortrows(sort_edges);
        f = find(diff([false;all(diff(B,1,1)==0,2);false])~=0);
        s = ones(length(f)/2,1);
        f1 = f(1:2:end-1); f2 = f(2:2:end);
        t = cumsum(accumarray([f1;f2+1],[s;-s],[size(B,1)+1,1]));
        sort_edges(ix(t(1:end-1)>0),:) = [];
        free_edges=sort_edges;
        if isempty(free_edges)
                watertight=1;
        else
                watertight=0;
        end
end