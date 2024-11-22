clc
clear all
close all

NormalVec = [h5read("../DFNGen.h5", "/NormalVec")]';
Location = [h5read("../DFNGen.h5", "/Location")]';
Phi = atan2(NormalVec(:, 2), NormalVec(:, 1));
Theta = acos(NormalVec(:, 3));

A = Phi < 0;
Phi(A) = Phi(A) + 2 * pi;

figure(1)
polarscatter(Phi, Theta, 's', 'filled'); rlim([0 0.5*pi]);
rticks([pi / 12, 2 * pi / 12, 3 * pi / 12, 4 * pi / 12, 5 * pi / 12, 6 * pi / 12 ]);
title(['Fractures', '''',' orientations']); hold on
set(gca,'rticklabel',[]);

% test normal vec
FractureVertices = [h5read("../DFNGen.h5", "/FractureVertices")]';
VerticesIndex = [h5read("../DFNGen.h5", "/VerticesIndex")]' + 1;

FractureRadius = zeros(size(NormalVec, 1), 1);
for i = 1:size(NormalVec, 1)
    Pnt1 = FractureVertices(VerticesIndex(i, 1), :);
    Pnt2 = FractureVertices(VerticesIndex(i, 1) + 1, :);
    % Pnt3 = FractureVertices(VerticesIndex(i, 1) + 2, :);

    S = Pnt2 - Pnt1;
    if dot(S, NormalVec(i, :)) > 1e-3
        error("Normal vector is not perpendicular to fracture")
    end

    FractureRadius(i) = norm(Location(i, :) - Pnt1);
end

% Size distribution
figure(2)
title("Fracture radius")
histogram(FractureRadius, 20);