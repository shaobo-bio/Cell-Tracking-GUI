function [Img Slice TimeZ Channel]= ULoadTiff_Z(filein)

% Function for loading in tif file from computer
% file in gives location of file for reading
% Reads in for 8 or 16 bit - for other formats saves as double

ImgProps    = imfinfo(filein);

Planes      = length(ImgProps);

Channel=0;
Slice=1;

if isfield (ImgProps,'ImageDescription')
    info=ImgProps(1).ImageDescription;
    
    if isempty(regexp(info,'slices'))
        Index=regexp(info,'channels')+9;
        Channel      = str2num(ImgProps(1).ImageDescription(Index));
    else
        Index=regexp(info,'slices')+7;
        Slice       = str2num(ImgProps(1).ImageDescription(Index:Index+1));
    end
    
end

if Channel == 0
    
    TimeZ       = Planes/Slice;
    
    ImgSize     = [ImgProps(1).Height, ImgProps(1).Width];
    
    % define data structure
    if(ImgProps(1).BitDepth == 8)
        Img     = uint8(zeros(ImgSize(1),ImgSize(2),Slice,TimeZ));
    elseif(ImgProps(1).BitDepth == 16)
        Img     = uint16(zeros(ImgSize(1),ImgSize(2),Slice,TimeZ));
    else
        disp('Warning, input image is neither 8 or 16 bit')
        Img     = zeros(ImgSize(1),ImgSize(2),Slice,TimeZ);
    end
    % size(Img)
    % load data
    for i = 1:TimeZ
        for j = 1: Slice
            
            Num=(i-1)*Slice+j;
            Img(:,:,j,i)  = imread(filein,Num);
            
        end
    end

    
else
    
    TimeZ       = Planes/Channel;
    
    ImgSize     = [ImgProps(1).Height, ImgProps(1).Width];
    
    % define data structure
    if(ImgProps(1).BitDepth == 8)
        Img     = uint8(zeros(ImgSize(1),ImgSize(2),TimeZ,Channel));
    elseif(ImgProps(1).BitDepth == 16)
        Img     = uint16(zeros(ImgSize(1),ImgSize(2),TimeZ,Channel));
    else
        disp('Warning, input image is neither 8 or 16 bit')
        Img     = zeros(ImgSize(1),ImgSize(2),TimeZ,Channel);
    end
    % size(Img)
    % load data
    for i = 1:TimeZ
        for j = 1: Channel
            
            Num=(i-1)*Channel+j;
            Img(:,:,i,j)  = imread(filein,Num);
            
        end
    end
end

Img=squeeze(Img);
