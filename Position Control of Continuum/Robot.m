classdef Robot
    % ROBOT Implements an n-link continuum robot
    
    properties
        nLinks % number of links
        OD     % nx1 vector contaning the OD of each link [m]
    end
    
    methods
        % Class constructor
        function self = Robot(nLinks, OD)
            self.nLinks = nLinks;
            self.OD = OD;
        end
         
        
        function [T, links] = fwkinematics(self, c)
            T = eye(4);
            links = cell(1,self.nLinks);
            
            for ii = 0 : self.nLinks - 1
               [Tii, link] = arckinematics(c(ii*3+1:ii*3+3));
               
               links{ii+1} = applytransform(link, T);
               T = T * Tii;
            end
        end
        
        
        function J = jacob(self, c)  
            

             jacobianSingleLink = @(k,phi,l) [cos(phi)*(cos(k*l)-1)/k^2 0 0;
                                              sin(phi)*(cos(k*l)-1)/k^2 0 0;
                                              -(sin(k*l)-k*l)/k^2 0 1;
                                              -l*sin(phi) 0 -k*sin(phi);
                                              l*cos(phi) 0 k*cos(phi);
                                              0 1 0];
             
             if nnz(jacobianSingleLink(c(1), c(2), c(3))) == 0
                 error('The Jacobian needs to be implemented.');
             end
             


             % add the jacobian for the first link
             J = zeros(6, 3 * self.nLinks);
             J(:,1:3) = jacobianSingleLink(c(1), c(2), c(3));
             
             % iteratively calculate and stack the jacobians for the other
             % links
             T = eye(4);
             
             for ii = 0 : self.nLinks - 2
                 % calculate the Adjoint transformation first
                 Tii = arckinematics(c(ii*3+1:ii*3+3));
                 T = T * Tii;
                 Ad = [T(1:3,1:3), skew(T(1:3,end)) * T(1:3,1:3);
                         zeros(3), T(1:3,1:3)];
                     
                 % then calculate the jacobian for the link
                 Jii = jacobianSingleLink(c(ii*3+4), c(ii*3+5), c(ii*3+6)); 
                 J(:,ii*3+4:ii*3+6) = Ad * Jii;
             end
        end
    end
end

