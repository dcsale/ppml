% compute and factorize the matrix for the bi-cubic interpolation in 2D 
% in ppm_gmm_kickoff.
%
% $Log: ppm_gmm_create_2dmatrix.m,v $
% Revision 1.1.1.1  2007/07/13 10:18:55  ivos
% CBL version of the PPM library
%
% Revision 1.2  2005/03/10 01:54:49  ivos
% Removed debug output.
%
% Revision 1.1  2005/03/10 01:53:00  ivos
% Initial check-in.
%
%-------------------------------------------------------------------------------
      A = zeros(16,16);
      for n=0:3,
	  nm1 = n-1;
	  % 0.0**n
	  ohn = 0.0;
	  if (n==0), ohn = 1.0; end;
	  % 0.0**(n-1)
	  ohnm1 = 0.0;
	  if (n==1), ohnm1 = 1.0; end;
	  for m=0:3,
	      mm1 = m-1;
	      % 0.0**m
	      ohm = 0.0;
	      if (m==0), ohm = 1.0; end;
	      % 0.0**(m-1)
	      ohmm1 = 0.0;
	      if (m==1), ohmm1 = 1.0; end;
                  % coefficient index
                  ind = 4*n+m+1;
                  % m --> x
                  % n --> y
                  %-------------------------------------------------------------
                  %  VALUES
                  %-------------------------------------------------------------
                  A(1,ind) = ohm*ohn;
                  A(2,ind) = ohn;
                  A(3,ind) = ohm;
                  A(4,ind) = 1.0;
                  %-------------------------------------------------------------
                  %  Dx
                  %-------------------------------------------------------------
                  A(5,ind) = m*ohmm1*ohn;
                  A(6,ind) = m*ohn;
                  A(7,ind) = m*ohmm1;
                  A(8,ind) = m;
                  %-------------------------------------------------------------
                  %  Dy
                  %-------------------------------------------------------------
                  A(9,ind) = n*ohm*ohnm1;
                  A(10,ind) = n*ohnm1;
                  A(11,ind) = n*ohm;
                  A(12,ind) = n;
                  %-------------------------------------------------------------
                  %  DxDy
                  %-------------------------------------------------------------
                  f         = m*n;
                  A(13,ind) = f*ohmm1*ohnm1;
                  A(14,ind) = f*ohnm1;
                  A(15,ind) = f*ohmm1;
                  A(16,ind) = f;
	      end;
	  end;

      [L,U,P] = lu(A);
      ALU = (L-eye(16))+U;
      ind = P*[1:1:16]';
      save 'gaga_ind.dat' ind -ASCII;
      ALUT = ALU';         % FORTRAN HAS COLUMN MAJOR
      save 'gaga_a.dat' ALUT -ASCII;

