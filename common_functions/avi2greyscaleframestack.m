function [ frame_stack, video ] = avi2greyscaleframestack( filename )
%avi2framestack Converts avi files to a 3D matrix
%
% Win's note:
% This function is no longer called from anywhere else in this repo.


video = VideoReader(filename);
video_data = video.read; % load the entire video into memory
frame_stack = zeros(video.Height, video.Width, video.NumberOfFrames,'uint8');
VideoFormat = video.VideoFormat;

switch VideoFormat
    case 'RGB24'
        for i=1:video.NumberOfFrames
            frame_stack(:,:,i) = rgb2gray(video_data(:,:,:,i));
        end
        
    case 'Grayscale'
        frame_stack = squeeze(video_data);
        
    otherwise
        warning('Unexpected VideoFormat. Modify the function accordingly.');
end

end

