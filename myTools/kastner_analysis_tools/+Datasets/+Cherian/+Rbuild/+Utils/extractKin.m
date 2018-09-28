function [ x ] = extractKin( mA, offset_Xpos, offset_Ypos )
% extracts kinematic data from madStruct
    
% mA : mad.analog ( Analog field of madStruct )    
% offset_Xpos : offset to put x position in terms of screen coords
% offset_Xpos : offset to put y position in terms of screen coords
% x : returns [ X Position ; Y Position ; X Velocity ; Y Velocity ]
    
% raw kinematic data
    tmp_Xpos = mA(19).waveform;
    tmp_Ypos = mA(20).waveform;
    
    % add offset to position
    Xposition = tmp_Xpos + offset_Xpos;
    Yposition = tmp_Ypos + offset_Ypos;
    
    % velocity data
    Xvelocity = mA(1).waveform;
    Yvelocity = mA(2).waveform;

    x = [Xposition Yposition Xvelocity Yvelocity];
end
