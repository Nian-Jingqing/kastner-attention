classdef Tensor
    properties
        t
    end

    methods
        % need to remove marginal means from dataset
        function out = removeMarginalMeans(t)
        % there is probably an efficient way to do this. I don't know it.
        % iterate over 3 dimensions, remove means
            for ix = 1:size(t.t, 1)
                tmean = mean( Tensor.vectorize( t.t(ix, :, :) ) );
                t.t(ix, :, :) = t.t(ix, :, :) - tmean;
            end
            
            for iy = 1:size(t.t, 2)
                tmean = mean( Tensor.vectorize( t.t( :, iy, :) ) );
                t.t( :, iy, :) = t.t( :, iy, :) - tmean;
            end
            
            for iz = 1:size(t.t, 3)
                tmean = mean( Tensor.vectorize( t.t(:, :, iz) ) );
                t.t( :, :, iz) = t.t( :, :, iz) - tmean;
            end
            
            out = t.t;
        end

        function out = calcCovarianceMatrices(t)
            tensor = t.t;
            % matrix names
            mnames = {'cx', 'cy', 'cz'};
            for ndim = 1:numel(mnames)
                numCounts = 0;
                out.(mnames{ndim}) = zeros( size(tensor, 1), size(tensor, 1) );
                for iy = 1:size(tensor, 2)
                    for iz = 1:size(tensor, 3)
                        out.(mnames{ndim}) = out.(mnames{ndim}) + ...
                            tensor(:, iy, iz) * tensor(:, iy, iz)';
                        numCounts = numCounts + 1;
                    end
                end
                out.(mnames{ndim}) = out.(mnames{ndim}) / numCounts;
                
                % blank out the diagonal elements?
                for nx = 1:size(tensor, 1)
                    out.(mnames{ndim})(nx, nx) = nan;
                end

                % shift dimensions for next calculation
                tensor = shiftdim( tensor, 1);
            end
        end
    end
end