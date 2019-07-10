function [ overwrite_flag ] = genTimeWarpDataFile( data, data_dir, data_name, overwrite_flag )
    data_file = [ data_name, '.mat' ];
    data_path = fullfile( data_dir, data_file );
    if exist( data_path, 'file' ) == 0 || overwrite_flag
        disp( 'INFO: Writing data file ...' )
        save( data_path, 'data' )
        overwrite_flag = false;
    else
        disp( 'INFO: Data file already created ...' )
        disp( 'INFO: Checking if input data matches loaded data ...' )
        S = load( data_path );
        x_check = isequal( S.data, data );
        
        if ~x_check
            disp( 'INFO: Input data did not match saved data. Overwriting data file ...' )
            save( data_path, 'data' )
        end
    end
end