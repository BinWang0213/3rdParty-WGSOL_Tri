function showEleSolution(nodes,u,v,edge_mid,u_edge,v_edge)
%% SHOWSOLUTION plots the velocity field on a triangular element in 2-D.
N=10;
[w1,w2] = meshgrid(linspace(0,1,N),linspace(0,1,N));
w3=1-w1-w2;

pts=[edge_mid];
U=[u_edge];
V=[v_edge];

%Create unifrom points within the triangular
for i=1:N
    for j=1:N
        %if(w1(i,j)<1e-8 | w2(i,j)<1e-8 | w3(i,j)<1e-8)
        if(w3(i,j)<0)
            continue;
        end
        w1(i,j)
        w2(i,j)
        w3(i,j)
        newpts=w1(i,j).*nodes(1,:) ...
              +w2(i,j).*nodes(2,:) ...
              +w3(i,j).*nodes(3,:);
        pts=[pts;newpts];

        newU=u(pts(end,1),pts(end,2));
        newV=v(pts(end,1),pts(end,2));
        U=[U;newU];
        V=[V;newV];
    end
end

%PlotBackground Mesh
figure(666);
h=patch(nodes(:,1),nodes(:,2),[0.5 0.9 0.45]);
h.EdgeColor='k';
h.FaceAlpha=0.5;
hold on;
%Plot interior velocity field
s=scatter(pts(:,1),pts(:,2),'filled'); %sample points
s.MarkerFaceColor='red';
quiver(pts(:,1),pts(:,2),U, V,'color','k');

%Plot edge velocity field
s=scatter(edge_mid(:,1),edge_mid(:,2),'filled'); %sample points
s.MarkerFaceColor='k';

view(2); axis equal; axis tight; axis off;
hold off;
end
% axis equal; axis tight;