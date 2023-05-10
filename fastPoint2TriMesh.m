function [distances,project_pt,outside]=fastPoint2TriMesh(inputs,pts,use_parallel)
    faces=inputs.faces;
    nodes=inputs.nodes;
    get_norm=0;
    if isfield(inputs,'face_mean_nodes')
        face_mean_nodes=inputs.face_mean_nodes;
    else
        get_norm=1;
    end
    if isfield(inputs,'face_normals');
        face_normals=inputs.face_normals;
    else
        get_norm=1;
    end
    
    get_tree=0;
    if isfield(inputs,'tree_model')
        tree_model=inputs.tree_model;
    else
        get_tree=1;
    end
    
    if get_norm==1
        [face_mean_nodes,face_normals]=getFaceCenterAndNormals(faces,nodes);
    end
    
    if get_tree==1
        tree_model=KDTreeSearcher(face_mean_nodes);
    end
    
    near_id=knnsearch(tree_model,pts);
    distances=dot(pts-face_mean_nodes(near_id,:),face_normals(near_id,:),2);
    
    signs=sign(distances);
    
    direction_vector=face_normals(near_id,:);
    project_pt=pts-(distances.*direction_vector);
    if use_parallel==1
        parfor count_pt=1:size(pts,1)
            dist=[0,0,0];
            q=zeros(3,3);
            project_pt_cur=project_pt(count_pt,:);


            node_ids=faces(near_id(count_pt),:);
            cor1=nodes(node_ids(1),:);
            cor2=nodes(node_ids(2),:);
            cor3=nodes(node_ids(3),:);
            cor4=project_pt_cur;
            
            check=f_check_inside_triangle( cor1,cor2,cor3,cor4);
            if check==0
                [dist(1),q(1,:)] = project_point_to_line_segment(cor1,cor2,cor4);
                [dist(2),q(2,:)] = project_point_to_line_segment(cor2,cor3,cor4);
                [dist(3),q(3,:)] = project_point_to_line_segment(cor1,cor3,cor4);
                [~,id]=min(dist);
                project_pt(count_pt,:)=q(id,:);
            end
        end

    else
        for count_pt=1:size(pts,1)
            project_pt_cur=project_pt(count_pt,:);


            node_ids=faces(near_id(count_pt),:);
            cor1=nodes(node_ids(1),:);
            cor2=nodes(node_ids(2),:);
            cor3=nodes(node_ids(3),:);
            cor4=project_pt_cur;
            
            check=f_check_inside_triangle( cor1,cor2,cor3,cor4);
            if check==0
                [dist(1),q(1,:)] = project_point_to_line_segment(cor1,cor2,cor4);
                [dist(2),q(2,:)] = project_point_to_line_segment(cor2,cor3,cor4);
                [dist(3),q(3,:)] = project_point_to_line_segment(cor1,cor3,cor4);
                [~,id]=min(dist);
                project_pt(count_pt,:)=q(id,:);
            end
        end
    end
    
    
    
    distances=vecnorm(pts-project_pt,2,2).*signs;
    
    outside=signs>=0;
end





%% helper functions

function [ores] = f_check_inside_triangle( cor1,cor2,cor3,cor4)
A=cor1;
B=cor2;
C=cor3;
P=cor4;   
    v0 = C - A;
    v1 = B - A;
    v2 = P - A;

    dot00 = dot(v0, v0);
    dot01 = dot(v0, v1);
    dot02 = dot(v0, v2);
    dot11 = dot(v1, v1);
    dot12 = dot(v1, v2);

    invDenom = 1 / (dot00 * dot11 - dot01 * dot01);
    u = (dot11 * dot02 - dot01 * dot12) * invDenom;
    v = (dot00 * dot12 - dot01 * dot02) * invDenom;
    ores =(u >= 0) && (v >= 0) && (u + v < 1);
end


function [dist,q] = project_point_to_line_segment(A,B,p)
      % returns q the closest point to p on the line segment from A to B 

      % vector from A to B
      AB = (B-A);
      % squared distance from A to B
      AB_squared = dot(AB,AB);
      if(AB_squared == 0)
        % A and B are the same point
        q = A;
      else
        % vector from A to p
        Ap = (p-A);
        % from http://stackoverflow.com/questions/849211/
        % Consider the line extending the segment, parameterized as A + t (B - A)
        % We find projection of point p onto the line. 
        % It falls where t = [(p-A) . (B-A)] / |B-A|^2
        t = dot(Ap,AB)/AB_squared;
        if (t < 0.0) 
          % "Before" A on the line, just return A
          q = A;
        else if (t > 1.0) 
          % "After" B on the line, just return B
          q = B;
        else
          % projection lines "inbetween" A and B on the line
          q = A + t * AB;
        end
        end
      end
      dist=norm(p-q);
end
