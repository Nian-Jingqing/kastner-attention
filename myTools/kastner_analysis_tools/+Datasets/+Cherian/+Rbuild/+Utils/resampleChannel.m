function [ resamp_data ] = resampleChannel( data, orig_sR, new_sR )
% Returns resampled data in a row matrix format
    flipFlag = false;
    if orig_sR < new_sR
        disp( 'WARNING: Upsampling data ...' )
    end
    
    if size(data,2) > size(data,1)
        data = data';
        flipFlag = true;
    end
    %[ b, a ] = butter( 2, new_sR / ( orig_sR / 2 ) );
    %filtered = filtfilt( b, a, data );
    %upsampleFactor = 2;
    %downsampleFactor = 5;
    %upsamp_data = upsample( filtered, upsampleFactor ); 
    %resamp_data = downsample( upsamp_data, downsampleFactor );
    resamp_data = resample( data, new_sR, orig_sR );
    %    if flipFlag
    %        resamp_data = resamp_data';
    %    end

end
