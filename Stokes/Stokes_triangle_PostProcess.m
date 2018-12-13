function errL2=Stokes_triangle_PostProcess(node,elem,h,info,pde,u)
% This is to solve the Stokes equation on triangular partition.
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

%Collect basic element info
NT = size(elem,1);
N  = size(node,1);
[elem2edge,edge] = dofedge(elem);
NE= size(edge,1);

%1. Get Edge flux and evaluate mass balance error
elem_edge_flux = zeros(NT,3);
for j=1:3
  elem_edge_flux(:,j) = u(elem2edge(:,j)).*info.Elem_edge_normal(:,j,1);
  elem_edge_flux(:,j) = elem_edge_flux(:,j)+ u(elem2edge(:,j)+NE).*info.Elem_edge_normal(:,j,2);
  elem_edge_flux(:,j) = elem_edge_flux(:,j).*info.Elem_edge_length(:,j);
end
elem_massbalance_error=sum(elem_edge_flux,2);


%2. Element-wise mass conservation velocity
err=[]; %Error term
u_centers=[];
u_exact_centers=[];
for ei=1:NT
    EleID=ei;
    EdgeIDs=elem2edge(EleID,:);

    %Node 13 18 17
    %Edge [17 18] [13 17] [13 18]
    %RT0 Edge-basis function
    Edge_u=u(elem2edge(EleID,:)); Edge_v=u(elem2edge(EleID,:)+NE);
    Edge_nx=info.Elem_edge_normal(EleID,:,1);
    Edge_ny=info.Elem_edge_normal(EleID,:,2);
    Edge_mid=info.Mid_Edge(EdgeIDs,1:2);

    RT0_Edge_u=Edge_u.*Edge_nx'+Edge_v.*Edge_ny';
    %Interpolate pts on elem center

    RT0_Basis=[];
    T=info.Elem_area(EleID);
    P=info.Elem_center(EleID,:); %Interpolate velocity at center
    for j=1:3
       ei=info.Elem_edge_length(EleID,j);
       Pi=node(elem(EleID,j),1:2);
       W=ei/2/T.*(P-Pi);
       RT0_Basis=[RT0_Basis,W'];
    end

    u_center_RT0=RT0_Basis*RT0_Edge_u;
    u_center_exact=pde.exactu(P);
    
    u_centers=[u_centers;u_center_RT0'];
    u_exact_centers=[u_exact_centers;u_center_exact];
    
    err=h*h.*[err;sum((u_center_RT0-u_center_exact').*(u_center_RT0-u_center_exact'))];
end
errL2=sqrt(sum(err));

%{
%Plot center uv
%Plot Element, debug
figure('rend','painters','pos',[10 10 800 600]);
set(gcf,'color','w');
showmesh(node,elem);
%findnode(node);
%findelem(node,elem);
axis on;
hold on;view(2);

quiver(info.Elem_center(:,1),info.Elem_center(:,2), ...
    u_centers(:,1),u_centers(:,2),'color','k');
quiver(info.Elem_center(:,1),info.Elem_center(:,2), ...
    u_exact_centers(:,1),u_exact_centers(:,2),'color','r');
legend('Mesh','RT0 Velocity','Exact Velocity',...
       'Orientation','horizontal','Location','northoutside')
hold off;
set(gca,'FontSize',15);
%}

end