function elasticksblktest()  
    tic;
    Li = 10;      % length of input zone  [cm]
    Lm = 80;      % length of middle zone [cm]
    Lo = 10;      % length of unloading   [cm]1
    
    Ks0 = 0.25;    % sieve conductivity [cm^2/(bar.s)]
    Kp = 2.5e-7;  % conductivity of membrane [cm/(bar.s)]
        
    r0 = 6e-4;     % radius of sieve [cm]
    A0 = pi*r0^2;   % cross-sectional area of sieve [cm]
    
    a  = 205;     % partial molar volume [cm^3/mol]
    
    psi = 0;      % water potential [bar]
        
    S = 4e-12;    % loading rate [mol/(cm.s)]
    ku = 5e-8;    % unloading rate [cm^2/s]
  
    
    T  = 297;     % temperature [K]
    R  = 83.141;  % [cm^3.bar/(mole.K)]
    
    %epsn=para(14);    %epsilon  Need to specify what this is, and what units we
    epsn=100; %right value in bar
    %         are using!!!

    dt=2;
    tprev=0;
    Nc=250;
    
    maxtime=2e5;
    
    % parameters related to blockage
    %
    t_block = maxtime;           % time blockage starts [s]
    t_open  = maxtime + 300;  % time blockage opens  [s]
    x_block = 45;             % spatial location of block [cm]
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    

    Lt = Li+Lm+Lo;         % total domain length
    dx = Lt/Nc;            % cell volume
    xc = dx*((1:Nc)'-0.5); % location of cell centers
    xe = [0; xc+dx/2];
 
    % make source and sink functions of space
    % S =  S*double( xc < Li ); %% PW constant source
    S = 2*S*(1-xc/Li).*double( xc < Li );  %modified source to smooth p
    ku = ku*double( xc > (Lt-Lo));
    
    % record index of blockage point
    %
    %k_block_m = max( find( xc < x_block ) );
    k_block_m = find( xc < x_block , 1, 'last' );
    k_block_p = k_block_m + 1;
    
    
    % matrix stuff  
    e = ones(Nc,1);
    B1 =spdiags([e e], -1:0,Nc,Nc);
    %B2 = spdiags([e 2*e e], -1:1,Nc,Nc);
    B3 = spdiags([e e], 0:1,Nc,Nc);
    
    
    
    % difference matrix for computing the velocity
    %   by differencing the pressure
    Dp = spdiags([-e e],[0 1],Nc-1,Nc);
    Dp = [zeros(1,Nc); Dp; zeros(1,Nc)];
    % modify the difference matrix when
    %   there is a blockage
    %
%     Dpblock = Dp;
%     Dpblock(k_block_p,:) = zeros(1,Nc);
    
    %difference matrix of volocity
    Du = spdiags([-ones(Nc+1,1) ones(Nc+1,1)],[0 1],Nc+1,Nc+1);
    
    
    % matrix for solving for pressure
    %
    
    
    % initalize
    %
    c = 0*xc;
    t = 0;
    r = r0*ones(Nc,1);     % radius of sieve [cm]
    A = A0*ones(Nc,1);   % cross-sectional area of sieve [cm]
    Ks = Ks0*ones(Nc,1);
    p0=0*ones(Nc,1);
    p=p0;
    u=zeros(Nc+1,1);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    recintv=5;
    recpoint=0:recintv:maxtime;
    k=1; %index of recorded time point
    numofrec=fix(maxtime/20);
    c_mat=zeros(Nc,numofrec);
    p_mat=zeros(Nc,numofrec);
    u_mat=zeros(Nc+1,numofrec);
    timearr=zeros(numofrec,1);
    xe_thin=zeros(6,1);
    xc_thin=zeros(6,1);
    
    cc=zeros(length(c),1000);
    pp=zeros(length(p),1000);
    uu=zeros(length(u),1000);
    AA=zeros(length(A),1000);
    tt=zeros(1000,1);
   
    
    
    
    irec=0;
    recflag=0;
    F(900) = struct('cdata',[],'colormap',[]);
    findex=1;
    while t<maxtime
        if( (t > t_block) && (t < t_open) )
            Ks(k_block_p-2)=Ks(k_block_p-2)*0.1;
            Ks(k_block_p-1)=Ks(k_block_p-1)*0.001;
            Ks(k_block_p)=Ks(k_block_p)*0.0001;
            Ks(k_block_p+1)=Ks(k_block_p+1)*0.001;
            Ks(k_block_p+2)=Ks(k_block_p+2)*0.1;
        end
        
        if( (t >= t_open) && (t < t_open+100) )
            Ks(k_block_p-1)=Ks(k_block_p-1)*(t-t_open)*0.01;
            Ks(k_block_p)=Ks(k_block_p)*(t-t_open)*0.01;
            Ks(k_block_p+1)=Ks(k_block_p+1)*(t-t_open)*0.01;
        end
        
        
        % choose a time step
        % flexible time step
        %dt = 0.8*dx/(max(abs(u)));

        
        % compute the upwind flux
        %
        cl = [0; c];
        cr = [c; 0];
        Al = [0; A];
        Ar = [A; 0];
        
        %J  = u.*cl.*double(u > 0) + u.*cr.*double( u < 0 );
        J  = u.*Al.*cl.*double(u > 0) + u.*Ar.*cr.*double( u < 0 );
        
        % solve for new concentration
        %
        %rhs = c - dt/dx*(J(2:end) - J(1:end-1)) + dt*S./A;
        %clast = c;
        %c  = rhs./(1+dt*ku./A);
        clast=c;
        CA=c.*A;
        Deltau=Du*u;
        Deltau=Deltau(1:end-1,:);
        rhs=S+CA.*(1/dt-Deltau/dx)-(J(2:end) - J(1:end-1))/dx;
        CA=rhs./((1/dt)+(ku./A));
        c=CA./A; %rediscover c
        
        
        % solve for a new pressure
        %
        rhs = p + dt.*(epsn./A).*(a.*(S-ku.*c)+2*pi*r*Kp.*(psi+R*T*c));

        % subdiagonal 
        %
        BA1 = 0.5*epsn*Ks/dx^2 .* (B1*A)./A; %Ks lifeng 10/23
        BA1(1:end-1)=BA1(2:end); %carefully figured out      
 
        % super diagonal
        %
        BA3=0.5*epsn*Ks/dx^2 .* (B3*A)./A; %Ks lifeng 10/23
        BA3(2:end)=BA3(1:end-1); %carefully figured out
           
        % main diagonal
        %
        %%%%BA2=B2*A*(-1/2).*((epsn*Kp)./(A.*dx^2));
        BA2=zeros(Nc,1);
        BA2 (2:end-1) = -BA1(1:end-2) - BA3(3:end);

        %%%%boundary condition here, bug fixed 11/09
        %BA2(1)=BA2(1)/2;
        BA2(1)=-epsn*Ks(1)/(A(1)*dx^2)*((A(2)+A(1))/2); 
        %BA2(end)=BA2(2)/2;
        BA2(end)=-epsn*Ks(end)/(A(end)*dx^2)*((A(end)+A(end-1))/2); 

           
        % form BA matrix from diagonals
        %
        %BA=spdiags([BA1 BA2 BA3],-1:1,Nc,Nc);
        
                 
            %BA=spdiags([BA1 BA2 BA3],-1:1,Nc,Nc);
            subdiag=-dt*BA1;
            diag=ones(Nc,1)+(epsn./A)*2*pi.*r*Kp*dt-BA2*dt;
            supdiag=-dt*BA3;
                   
            D  = -1/dx * Dp;       

        
        
        %%% RDG 10/05/12
        %
        % it looks to me like BA has some positive eigenvalues;  I don't
        % think this is right, and we need to check this carefully!!!!
        %
        % just glacing at a corner of the matrix, I can see that there is
        %   something wrong at the endpoints it looks like the first and
        %   last row are wrong
        %
        % Perhaps we fixed this in the afternoon.  Lifeng will check 
        %  it carefully
        %
        
        %test =spdiags(1+(epsn./A)*2*pi.*r*Kp*dt,0,Nc,Nc) -dt.*BA;
        lhs=spdiags([subdiag diag supdiag],-1:1,Nc,Nc);
        p = lhs\rhs;
    
        %%plot(xc,p,'ro');
        %%title(sprintf('time = %g',t));
        %%pause(0.1);
        %%%        keyboard
    
        
        Kstmp=B1*Ks*0.5;
        Kstmp=[0;Kstmp(1:end-1);0];
        u  = D*p.*Kstmp;
        A = A0*ones(Nc,1).*exp((p-p0)/epsn);
        r = sqrt(A/pi);
        Ks=(r/r0).^4*Ks0;
        t = t + dt;
        % update the time
        
        %%adjust time step
        if dt>0.1*dx/(max(abs(u)))
            %dt=0.1*double(0.1*dx/(max(abs(u)))>0.1)+0.01*double(0.1*dx/(max(abs(u)))<0.1);
            dt=0.05*dx/(max(abs(u)));
        elseif dt<0.05*dx/(max(abs(u))) && dt<2
            %dt=1*double(0.1*dx/(max(abs(u)))>1)+0.1*double(0.1*dx/(max(abs(u)))<1);
            dt=min(0.07*dx/(max(abs(u))),2);
        end
        if dt>0.8*dx/(max(abs(u)))
            checkstep=dt/(dx/(max(abs(u))));      
            warning('checkstep=%f dt=%d',checkstep,dt);
        end
        %dt=0.1;
        
        if ( t>=recpoint(k)) %%%need adapatable code
            k=k+1;
            %disp(t);
            if recflag==1
                irec=irec+1;
                timearr(irec)=t;
                c_mat(:,irec)=c;
                p_mat(:,irec)=p;
                u_mat(:,irec)=u;
            end
        end
        
        if( max( abs( 1 - c./clast )/dt) < 1e-5)
            if recflag==0               
                blocktime=t
                recflag=1;
                t_block=t+600;
                t_open=t_block+600;
            elseif t>t_open+420
                sprintf('blockage happened at %d',t_block)
                sprintf('2nd steady state at %d, program ended',t)
                break;
            end
        end


        %if ( mod(t,5)<0.01 && t>1.5e+04)
        %if ( mod(t,5)<0.01) 
        if ( t-tprev>2 && t>t_block-60 && t<t_open+600)
            tprev=t;

            
            hf=figure(1);
            %set(hf,'Position',[0 0 1400 700])
            subplot(4,1,1);
            plot(xc,c,'r.-');
            ylim([0 9e-4])
            set(gca,'fontsize',12);
            title((strcat('Time since blockage=',num2str(t-t_block))));
            ylabel('concentration');

            subplot(4,1,2);
            plot(xc,p,'b.-');
            ylim([0 25])
            set(gca,'fontsize',12);
            ylabel('pressure');

            subplot(4,1,3);
            plot(xe,60*u,'k.-');
            xlim([0 xe(end)])
            %ylim([0 4])
            set(gca,'fontsize',12);
            ylabel('velocity');
            
            subplot(4,1,4);
            plot(xc,A,'c.-');
            ylim([0 1.7e-6])
            set(gca,'fontsize',12);
            ylabel('sectional area');
            F(findex) = getframe(hf);
            
            pause(0.1);
            tt(findex)=t-t_block;
            cc(:,findex)=c;
            pp(:,findex)=p;
            uu(:,findex)=u;
            AA(:,findex)=A;
            findex=findex+1;
       end
        
    end
    
    
    toc;
    if recflag==0
        error('blockage has not happened, please increase maxtime');
    end

    %%%%% animation stuff
%     field1 = 'time';  value1 = tt;
%     field2 = 'c';  value2 = cc;
%     field3 = 'p';  value3 = pp;
%     field4 = 'u';  value4 = uu;
%     field5 = 'A';  value5 = AA;
%     s = struct(field1,value1,field2,value2,field3,value3,field4,value4,field5,value5,'xc',xc,'xe',xe);
%     return
%%%%%%%%%%%%%%%%%%%%%

    timearr(irec:end)=[];
    c_mat(:,irec:end)=[];
    p_mat(:,irec:end)=[];
    u_mat(:,irec:end)=[];
%     [file,path] = uiputfile('*.xls','Save outpu As excel');
%     output=strcat(path,file);
    ithin=0;
    numofrec=length(c_mat);
    c_mat_thin=zeros(6,numofrec);
    p_mat_thin=zeros(6,numofrec);
    u_mat_thin=zeros(6,numofrec);
    for irow=[38 63 88 138 163 188]
        ithin=ithin+1;
        xc_thin(ithin,:)=xc(irow,:);
        xe_thin(ithin,:)=xe(irow,:);
        c_mat_thin(ithin,:)=c_mat(irow,:);
        p_mat_thin(ithin,:)=p_mat(irow,:);
        u_mat_thin(ithin,:)=u_mat(irow,:);
    end
    keyboard
    
%     xlswrite(output,xc_thin',1,'B1');
%     xlswrite(output,timearr,1,'A2');
%     xlswrite(output,c_mat_thin',1,'B2');
%     
%     xlswrite(output,xc_thin',2,'B1');
%     xlswrite(output,timearr,2,'A2');
%     xlswrite(output,p_mat_thin',2,'B2');
%     
%     xlswrite(output,xe_thin',3,'B1');
%     xlswrite(output,timearr,3,'A2');
%     xlswrite(output,u_mat_thin',3,'B2');
    
%     xlswrite(output,xc',4,'B1');
%     xlswrite(output,timearr,4,'A2');
%     xlswrite(output,c_mat',4,'B2');
%     
%     xlswrite(output,xc',5,'B1');
%     xlswrite(output,timearr,5,'A2');
%     xlswrite(output,p_mat',5,'B2');
%     
%     xlswrite(output,xe',6,'B1');
%     xlswrite(output,timearr,6,'A2');
%     xlswrite(output,u_mat',6,'B2');

end