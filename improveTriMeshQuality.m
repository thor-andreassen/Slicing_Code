function [nodes,faces]=improveTriMeshQuality(faces,vertices,num_neighbors,max_iters,tol)
counter=1;
conv=Inf;
nodes=vertices;
for count_node=1:size(nodes,1)
        [near_nodes_id{count_node},ring_nodes_layered{count_node}]=findKRingNodes(faces,count_node,num_neighbors);
end

while counter<=max_iters && conv>tol
    current_order=1:size(nodes,1);
    current_order = current_order(randperm(length(current_order)));
    conv=-Inf;
    count_inner=1;
    for count_node=current_order
        %                [new_node_loc,disp_err]=shiftNodeLocation(faces,nodes,count_node,num_neighbors);
        try
            [new_node_loc,disp_err]=shiftNodeLocation(faces,nodes,count_node,...
                near_nodes_id{count_node},ring_nodes_layered{count_node});
        catch
            new_node_loc=nodes(count_node,:);
        end
        
        if disp_err>=conv
            conv=disp_err;
        end
        nodes(count_node,:)=new_node_loc;
        count_inner=count_inner+1;
        if mod(count_inner,100)==0
            disp(num2str(count_inner));
        end
    end
    counter=counter+1;
end