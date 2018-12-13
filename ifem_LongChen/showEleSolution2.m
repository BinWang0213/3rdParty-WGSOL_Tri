function showEleSolution2(nodes,uv,edge_mid,u_edge,v_edge)
    %% SHOWSOLUTION plots the velocity field on a triangular element in 2-D.
    N=10;
    [w1,w2] = meshgrid(linspace(0,1,N),linspace(0,1,N));
    w3=1-w1-w2;
    
    pts=[];
    U=[];
    V=[];
    
    %Create unifrom points within the triangular
    for i=1:N
        for j=1:N
            %if(w1(i,j)<1e-8 | w2(i,j)<1e-8 | w3(i,j)<1e-8)
            if(w3(i,j)<0)   
                continue;
            end
            newpts=w1(i,j).*nodes(1,:) ...
                  +w2(i,j).*nodes(2,:) ...
                  +w3(i,j).*nodes(3,:);
            pts=[pts;newpts];
            newUV=uv(pts(end,1),pts(end,2));
            if(abs(newUV(1))<1e-10 && abs(newUV(2))<1e-10)
                newUV=[0 0];
            end
            U=[U;newUV(1)];
            V=[V;newUV(2)];
        end
    end
    
    %PlotBackground Mesh
    figure(777);
    h=patch(nodes(:,1),nodes(:,2),[0.5 0.9 0.45]);
    h.EdgeColor='k';
    h.FaceAlpha=0.5;
    hold on;
    %Plot interior velocity field
    s=scatter(pts(:,1),pts(:,2),'filled'); %sample points
    s.MarkerFaceColor='red';
    quiver(pts(:,1),pts(:,2),U, V,'color','red');
    
    %Plot edge velocity field
    pts=[edge_mid];
    U=[u_edge];
    V=[v_edge];
    quiver(pts(:,1),pts(:,2),U, V,'color','k');
    s=scatter(edge_mid(:,1),edge_mid(:,2),'filled'); %sample points
    s.MarkerFaceColor='k';
    
    view(2); axis equal; axis tight; axis on;
    hold off;
    end
    % axis equal; axis tight;