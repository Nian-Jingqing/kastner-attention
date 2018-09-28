function ordering = getCrossCorrOrdering( matrix, cutoff )

if ~exist('cutoff', 'var')
    error('getCrossCorrOrdering: need to specify a cutoff percentile');
end

% set diagonals to nan
for id = 1:size( matrix, 1 )
    matrix( id, id ) = nan;
end


matrix2 = matrix;
%  make the matrix sparse
matrix2( matrix  <= ...
        prctile( matrix(:), cutoff ) ) = 0;
matrix2( matrix2 > 0) = 1;



% set diagonals to 0
for id = 1:size( matrix, 1 )
    matrix2( id, id ) = 0;
end



% use maximal clique alg

orderingTmp = Utils.maximalCliques( matrix2 );


%% go through the first N, cliques, pull out the ids, stack them
list = [];
for nc = 1:size( orderingTmp, 2 )
    list = vertcat(list(:), find( orderingTmp( :, nc ) ) );
end
[c, ia, ic ] = unique( list, 'first' );
[~, inds] = sort(ia);
ordering = inds;
