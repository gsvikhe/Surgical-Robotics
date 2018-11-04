classdef Wrist

   properties
      ID                % Inner Diameter in mm
      OD                % Outer diameter in mm
      n_notch         % # of notches
      notch_struct      % = struct('h', 'x', 'w', 'theta');  
      
      %{ 
      notch_struct is a structure with the following fields
      h - height of notch
      x - distance between notches                                  
      w - width of notches             
      orient - orientation of notch
      %}
   end
   
   methods       
      function obj = Wrist(ID, OD, n_notch, notch_struct)
               obj.ID = ID;
               obj.OD = OD;
               obj.n_notch = n_notch;
               obj.notch_struct = notch_struct;
      end
      
      function trans_mat = FwKin(obj, configuration)
            % Read wrist design
            n  = obj.n_notch;                                              % total number of notch
            ro = obj.OD/2 ;                                                % outer radius of tube 
            ri = obj.ID/2 ;                                                % inner radius of tube 

            x = obj.notch_struct.x ;                                       % Length between notch 
            h = obj.notch_struct.h ;                                       % height of the notch 
            w = obj.notch_struct.w ;                                       % notch depth 
            notch_orient = [obj.notch_struct.orient(1)] * pi / 180;        % notch orientation in radians.
            
            
            % Get the wrist configuration
            tendon_disp = configuration(1) ;                   % displacement
            wrist_rot   = configuration(2) * pi / 180;         % wrist rotation in radians
            wrist_adv   = configuration(3) ;                   % Wrist advancement in the Z-axis 

            temp = w-ro; % temp variable
            % Transformation Matrix corresponding to rotation about the z-axis
            T_rot_z = @(alpha) [cos(alpha) -sin(alpha) 0 0; 
                                sin(alpha) cos(alpha)  0 0; 
                                0          0           1 0; 
                                0          0           0 1];
 
            % Transformation Matrix corresponding to translation along z-axis
            T_trans_z = @(a) [1 0 0 0;
                              0 1 0 0;
                              0 0 1 a;
                              0 0 0 1];     
                          
                          
            % Transformation Matrix corresponding to rotation
            % about the x-axis and translation in y and z. From York Paper. Equation 6
            T_notch = @(kappa,s) [1 0 0 0;
                                   0 cos(kappa*s) -sin(kappa*s) (cos(kappa*s)-1)/kappa;
                                   0 sin(kappa*s) cos(kappa*s) sin(kappa*s)/kappa;
                                   0 0 0 1];    
                               

            % Calculating Ybar from York Paper. Equation 1. 
            phi_outer = 2 * acos((temp)/ro);
            phi_inner = 2 * acos((temp)/ri);
            
            A_outer = ro^2 * (phi_outer-sin(phi_outer))/2;
            A_inner = ri^2 * (phi_inner-sin(phi_inner))/2;
            
            ybar_outer = 4 * ro * sin(0.5 * phi_outer)^3/ 3 *( phi_outer - sin(phi_outer));
            ybar_inner = 4 * ro * sin(0.5 * phi_inner)^3/ 3 * (phi_inner - sin(phi_inner));
            
            % Using the small angle approximation to calculate ybar
            ybar = (ybar_outer * A_outer - ybar_inner * A_inner)/(A_outer- A_inner)

            % Kappa of Single Cutout based on small angle approximation.
            % York Paper Equation 4
            kappa = tendon_disp/(h * (ri + ybar)- tendon_disp * ybar);

            % Arc length of single cutout
            s = h / ( 1 + ybar * kappa);

            % Initialize the transformation matrix
            T = repmat(eye(4), 1, 1, n + 2);                    % n+2 because number of notches + base + tip

            % Calculate Transformation Matrix based on initial advancement
            % and rotation
            T(:,:,2) = T_trans_z(wrist_adv) * T_rot_z(wrist_rot); 
            
            % Iterate on the cutouts and calculate the transformations at each cutout
            for i = 3 : (n + 2)
                T(:,:,i) = T(:,:,i-1) * T_rot_z(notch_orient) * T_notch(kappa, s) * T_trans_z(x);    
            end
            % Extract the tube pose and return
            trans_mat = T;
            theta_max = h/(ro+ybar)*180/pi
      end
      
      
      function max_strain = calcmaxstrain(obj, configuration)
          % Read wrist design
            n  = obj.n_notch;                                % total number of notches
            ro = obj.OD/2 ;                           % outer radius of tube 
            ri = obj.ID/2 ;                           % inner radius of tube 

            x = obj.notch_struct(1).x ;                      % Length between cuts 
            h = obj.notch_struct(1).h ;                      % height of the cutouts 
            w = obj.notch_struct(1).w ;                       % cut depth 
            notch_orient = [obj.notch_struct.orient] * pi / 180;     % notch orientation in radians.
            temp = w-ro; 
            
            % Get the wrist configuration
            tendon_disp = configuration(1) ;            % displacement 
            wrist_rot   = configuration(2) ;            % wrist rotation
            wrist_adv   = configuration(3) ;            % Wrist advancement in the Z-axis 
            
            % Calculating Ybar from York Paper. Equation 1. 
            phi_outer = 2 * acos((temp)/ro);
            phi_inner = 2 * acos((temp)/ri);
            
            A_outer = ro^2 * (phi_outer-sin(phi_outer))/2;
            A_inner = ri^2 * (phi_inner-sin(phi_inner))/2;
            
            ybar_outer = (4 * ro * sin(0.5 * phi_outer)^3)/ (3 * ( phi_outer - sin(phi_outer)));
            ybar_inner = (4 * ro * sin(0.5 * phi_inner)^3)/ (3 * (phi_inner - sin(phi_inner)));
            
            % Using the small angle approximation to calculate ybar
            ybar = (ybar_outer * A_outer - ybar_inner * A_inner)/(A_outer- A_inner);

            % Kappa of Single Cutout based on small angle approximation.
            % York Paper Equation 4
            kappa = tendon_disp/(h * (ri + ybar)- tendon_disp * ybar);
            
            max_strain = (kappa*(ro-ybar)/(1+ybar*kappa));
            if max_strain > 0.08 | max_strain<-0.08
                f = msgbox(strcat('The max strain exceeds the recoverable strain of Nitinol at tendon_disp equal to','  ',num2str(tendon_disp)), 'Warning','warn');
            end
      end
      
      function trans_mat = FwKin2(obj,configuration)
          if length(configuration)==4
            % Read wrist design
            n  = obj.n_notch;                                             % total number of notches
            ro = obj.OD/2;                                                % outer radius of tube 
            ri = obj.ID/2;                                                % inner radius of tube 

            x = obj.notch_struct(1).x ;                                   % Length between notches 
            h = obj.notch_struct(1).h ;                                   % height of the notches 
            w = obj.notch_struct(1).w ;                                   % notch depth 
            
            
            notch_orient_1 = [obj.notch_struct.orient(1)] * pi / 180;     % notch orientation in radians.
            notch_orient_2 = [obj.notch_struct.orient(2)] * pi / 180;
            
            % Get the wrist configuration
            tendon_disp_1 = configuration(1); 
            tendon_disp_2 = configuration(2);  
            wrist_rot   = configuration(3) * pi / 180;
            wrist_adv   = configuration(4); 

            temp = w-ro; % temp variable 

            
            % Transformation Matrix corresponding to rotation about the x-axis 
            % and translation in y and z. From York Paper. Equation 6
            T_notch = @(kappa,s) [1 0 0 0;
                                   0 cos(kappa*s) -sin(kappa*s) (cos(kappa*s)-1)/kappa;
                                   0 sin(kappa*s) cos(kappa*s) sin(kappa*s)/kappa;
                                   0 0 0 1];       

            % Transformation Matrix corresponding to rotation about the z-axis
            T_rot_z = @(alpha) [cos(alpha) -sin(alpha) 0 0; 
                                sin(alpha) cos(alpha)  0 0; 
                                0          0           1 0; 
                                0          0           0 1];
            % Transformation Matrix corresponding to translation along z-axis
            T_trans_z = @(a) [1 0 0 0;
                              0 1 0 0;
                              0 0 1 a;
                              0 0 0 1];     

            % Calculating Ybar from York Paper. Equation 1. 
            phi_outer = 2 * acos((temp)/ro);
            phi_inner = 2 * acos((temp)/ri);
            
            A_outer = ro^2 * (phi_outer-sin(phi_outer))/2;
            A_inner = ri^2 * (phi_inner-sin(phi_inner))/2;
            
            ybar_outer = 4 * ro * sin(0.5 * phi_outer)^3/ 3 * (phi_outer - sin(phi_outer));
            ybar_inner = 4 * ro * sin(0.5 * phi_inner)^3/ 3 * (phi_inner - sin(phi_inner));
            
            ybar = (ybar_outer * A_outer - ybar_inner * A_inner)/(A_outer- A_inner);
            
            % Kappa of Single Cutout based on small angle approximation.
            % York Paper Equation 4
            kappa1 = tendon_disp_1/(h * (ri + ybar)- tendon_disp_1 * ybar);
            kappa2 = tendon_disp_2/(h * (ri + ybar)- tendon_disp_2 * ybar);

            % Get arc length of single cutout
            s_1 = h / ( 1 + ybar * kappa1);
            s_2 = h / ( 1 + ybar * kappa2);

            % Initialize transformation matrix
            T = repmat(eye(4), 1, 1, n + 2);

            % Calculate transformation matrix based on initial advancement and rotation
            T(:,:,2) = T_trans_z(wrist_adv) * T_rot_z(wrist_rot); 
            
            % Calculate the transformations at each cutout
            for i = 3 : n
                T(:,:,i) = T(:,:,i-1) * T_rot_z(notch_orient_1) * T_notch(kappa1, s_1) * T_trans_z(x);    
            end
            for i = n+1 : n+2
                T(:,:,i) = T(:,:,i-1) * T_rot_z(notch_orient_2) * T_notch(kappa2, s_2) * T_trans_z(x);    
            end
            
%             T(:,:,3) = T(:,:,2) * T_rot_z(notch_orient_1) * T_notch(kappa1, s_1) * T_trans_z(x);    
%             T(:,:,4) = T(:,:,3) * T_rot_z(notch_orient_2) * T_notch(kappa2, s_2) * T_trans_z(x); 
%             T(:,:,5) = T(:,:,4) * T_rot_z(notch_orient_1) * T_notch(kappa1, s_1) * T_trans_z(x);
%             T(:,:,6) = T(:,:,5) * T_rot_z(notch_orient_2) * T_notch(kappa1, s_1) * T_trans_z(x); 
            
            % The Final Transformation Matrix is thus
            trans_mat = T;
          else
              disp("Put in 4 values in configuration array. The first two representing the tendon displacements")
          end
      end
      
      
      
      function trans_mat = FwKin3(obj,configuration)
          if length(configuration)==6
            % Read wrist design
            n  = obj.n_notch;                                             % total number of notches
            ro = obj.OD/2;                                                % outer radius of tube 
            ri = obj.ID/2;                                                % inner radius of tube 

            x = obj.notch_struct(1).x ;                                   % Length between notches 
            h = obj.notch_struct(1).h ;                                   % height of the notches 
            w = obj.notch_struct(1).w ;                                   % notch depth 
            
            
            notch_orient_1 = [obj.notch_struct.orient(1)] * pi / 180;     % notch orientation in radians.
            notch_orient_2 = [obj.notch_struct.orient(2)] * pi / 180;
            notch_orient_3 = [obj.notch_struct.orient(3)] * pi / 180;
            notch_orient_4 = [obj.notch_struct.orient(4)] * pi / 180;
            % Get the wrist configuration
            tendon_disp_1 = configuration(1); 
            tendon_disp_2 = configuration(2);  
            tendon_disp_3 = configuration(3);
            tendon_disp_4 = configuration(4);
            
            wrist_rot   = configuration(3) * pi / 180;
            wrist_adv   = configuration(4); 

            temp = w-ro; % temp variable 

            
            % Transformation Matrix corresponding to rotation about the x-axis 
            % and translation in y and z. From York Paper. Equation 6
            T_notch = @(kappa,s) [1 0 0 0;
                                   0 cos(kappa*s) -sin(kappa*s) (cos(kappa*s)-1)/kappa;
                                   0 sin(kappa*s) cos(kappa*s) sin(kappa*s)/kappa;
                                   0 0 0 1];       

            % Transformation Matrix corresponding to rotation about the z-axis
            T_rot_z = @(alpha) [cos(alpha) -sin(alpha) 0 0; 
                                sin(alpha) cos(alpha)  0 0; 
                                0          0           1 0; 
                                0          0           0 1];
            % Transformation Matrix corresponding to translation along z-axis
            T_trans_z = @(a) [1 0 0 0;
                              0 1 0 0;
                              0 0 1 a;
                              0 0 0 1];     

            % Calculating Ybar from York Paper. Equation 1. 
            phi_outer = 2 * acos((temp)/ro);
            phi_inner = 2 * acos((temp)/ri);
            
            A_outer = ro^2 * (phi_outer-sin(phi_outer))/2;
            A_inner = ri^2 * (phi_inner-sin(phi_inner))/2;
            
            ybar_outer = 4 * ro * sin(0.5 * phi_outer)^3/ 3 * (phi_outer - sin(phi_outer));
            ybar_inner = 4 * ro * sin(0.5 * phi_inner)^3/ 3 * (phi_inner - sin(phi_inner));
            
            ybar = (ybar_outer * A_outer - ybar_inner * A_inner)/(A_outer- A_inner);
            
            % Kappa of Single Cutout based on small angle approximation.
            % York Paper Equation 4
            kappa1 = tendon_disp_1/(h * (ri + ybar)- tendon_disp_1 * ybar);
            kappa2 = tendon_disp_2/(h * (ri + ybar)- tendon_disp_2 * ybar);
            kappa3 = tendon_disp_3/(h * (ri + ybar)- tendon_disp_3 * ybar);
            kappa4 = tendon_disp_4/(h * (ri + ybar)- tendon_disp_4 * ybar);

            
            % Get arc length of single cutout
            s_1 = h / ( 1 + ybar * kappa1);
            s_2 = h / ( 1 + ybar * kappa2);
            s_3 = h / ( 1 + ybar * kappa3);
            s_4 = h / ( 1 + ybar * kappa4);

            % Initialize transformation matrix
            T = repmat(eye(4), 1, 1, n + 2);

            % Calculate transformation matrix based on initial advancement and rotation
            T(:,:,2) = T_trans_z(wrist_adv) * T_rot_z(wrist_rot); 
            
            % Calculate the transformation matrix
            T(:,:,3) = T(:,:,2) * T_rot_z(notch_orient_1) * T_notch(kappa1, s_1) * T_trans_z(x);    
            T(:,:,4) = T(:,:,3) * T_rot_z(notch_orient_2) * T_notch(kappa2, s_2) * T_trans_z(x);
            T(:,:,5) = T(:,:,4) * T_rot_z(notch_orient_3) * T_notch(kappa3, s_3) * T_trans_z(x);
            T(:,:,6) = T(:,:,5) * T_rot_z(notch_orient_4) * T_notch(kappa4, s_4) * T_trans_z(x);
            T(:,:,7) = T(:,:,6) * T_rot_z(notch_orient_1) * T_notch(kappa1, s_1) * T_trans_z(x);    
            T(:,:,8) = T(:,:,7) * T_rot_z(notch_orient_2) * T_notch(kappa2, s_2) * T_trans_z(x);
            T(:,:,9) = T(:,:,8) * T_rot_z(notch_orient_3) * T_notch(kappa3, s_3) * T_trans_z(x);
            T(:,:,10) = T(:,:,9) * T_rot_z(notch_orient_4) * T_notch(kappa4, s_4) * T_trans_z(x);
            
            % The Final Transformation Matrix is thus
            trans_mat = T;
          else
              disp("Put in 6 values in configuration array. The first two representing the tendon displacements")
          end
      end
   end    
end