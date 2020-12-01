function [MOVE_PointCloud] =...
    Create_MOVEPointCloud(x_max, x_min, z_max, z_min, resolution)
%   INPUT:
%   x_max, x_min, z_max, z_min: integers of bounding box in [km]
%   resolution: MOVE point cloud resolution in [km]

%   OUTPUT:
%   columns: x,y,z, color ID, unique ID
%
%   Paul R. Eizenh?fer, PhD
%   University of Pittsburgh, peizen@pitt.edu

[X,Z] = meshgrid(x_min:resolution:x_max, z_min:resolution:z_max);

MOVE_PointCloud(:,1) = X(:);
MOVE_PointCloud(:,2) = 0;
MOVE_PointCloud(:,3) = Z(:);
MOVE_PointCloud(:,4) = 0; 
for i = 1:length(MOVE_PointCloud(:,3))
    MOVE_PointCloud(i, 5) = i;
end