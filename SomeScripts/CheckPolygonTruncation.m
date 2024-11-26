clc
clear all
close all

stlFile = '../Box3D.stl';

% Read the STL file
model = stlread(stlFile);

% Extract faces and vertices
faces = model.ConnectivityList;
vertices = model.Points;

% Display the STL model
figure(1);
patch('Vertices', vertices, 'Faces', faces, ...
      'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0);


OriginalP = [0.100408, 0.682157, 0.361363;  0.30806, 0.825049, 0.176795;  0.0887146, 0.992077, 0.323736;];

StructP = [1,2,3];

Pkk = [0.0887146, 0.992077, 0.323736;  0.100408, 0.682157, 0.361363;  0.30806, 0.825049, 0.176795;  0.30806, 0.825049, 0.176795;  0.0887146, 0.992077, 0.323736;  0.100408, 0.682157, 0.361363];

figure(1)
view(3)
patch('Vertices', OriginalP, 'Faces', StructP, 'FaceVertexCData', zeros(size(StructP, 1), 1), ...
    'FaceColor', 'flat', 'EdgeAlpha', 1, 'facealpha', 0); hold on

for i = 1:size(Pkk, 1)
    scatter3(Pkk(i, 1), Pkk(i, 2), Pkk(i, 3), '*'); hold on
end