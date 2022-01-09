function [new_node_loc,disp]=shiftNodeLocation(faces,nodes,num_node,...
        near_nodes_id,ring_nodes_layered)
        near_nodes=nodes(near_nodes_id,:);
%         [near_nodes_id,ring_nodes_layered]=findKRingNodes(faces,num_node,num_neighbors);

        vecz=findVertexNormalFromMesh(faces,nodes,num_node);
        temp_vecx=near_nodes(2,:)-near_nodes(1,:);
        temp_vecx=temp_vecx/norm(temp_vecx);
        vecy=cross(vecz,temp_vecx)/norm(cross(vecz,temp_vecx));
        vecx=cross(vecy,vecz)/norm(cross(vecy,vecz));

        rotmat=[vecx', vecy', vecz'];
        transvec=[nodes(num_node,:)'];

        LFCS_to_GCS_TransMat=eye(4);
        LFCS_to_GCS_TransMat(1:3,1:3)=rotmat;
        LFCS_to_GCS_TransMat(1:3,4)=transvec;
        GCS_to_LFCS_TransMat=inv(LFCS_to_GCS_TransMat);

        near_nodes_surf_project=GCS_to_LFCS_TransMat*[near_nodes';ones(1,size(near_nodes,1))];
        near_nodes_surf_project=near_nodes_surf_project(1:3,:)';


        nearest_ring_id=ring_nodes_layered{2};
        nearest_ring_nodes=nodes(nearest_ring_id,:);
        near_nodes_polygon_project=GCS_to_LFCS_TransMat*[nearest_ring_nodes';ones(1,size(nearest_ring_nodes,1))];
        near_nodes_polygon_project=near_nodes_polygon_project(1:3,:)';

        temp_res=max_inscribed_circle(near_nodes_polygon_project(:,1),near_nodes_polygon_project(:,2));
        new_node_loc_proj(1)=temp_res.xc;
        new_node_loc_proj(2)=temp_res.yc;
        [quadric_coefs,~]=fitQuadricSurf(near_nodes_surf_project(:,1),near_nodes_surf_project(:,2),near_nodes_surf_project(:,3));

        quadric_vec=[new_node_loc_proj(1).^2,...
                new_node_loc_proj(1).*new_node_loc_proj(2),...
                new_node_loc_proj(2).^2,...
                new_node_loc_proj(1),...
                new_node_loc_proj(2),...
                1];
        new_node_loc_proj(3)=quadric_vec*quadric_coefs;
        new_node_loc_GCS=LFCS_to_GCS_TransMat*[new_node_loc_proj';ones(1,size(new_node_loc_proj,1))];
        new_node_loc_GCS=new_node_loc_GCS(1:3)';
        new_node_loc=new_node_loc_GCS;
        disp=norm(new_node_loc-near_nodes(1,:));