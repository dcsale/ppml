% compute and factorize the matrix for the tri-cubic interpolation in 3D 
% in ppm_gmm_kickoff.
%
% $Log: ppm_gmm_create_3dmatrix.m,v $
% Revision 1.1.1.1  2007/07/13 10:18:55  ivos
% CBL version of the PPM library
%
% Revision 1.1  2005/03/10 01:53:00  ivos
% Initial check-in.
%
%-------------------------------------------------------------------------------
      A = zeros(64,64);
      for p=0:3,
          pm1 = p-1;
          % 0.0**p
          ohp = 0.0;
	  if (p==0), ohp = 1.0; end;
          % 0.0**(p-1)
          ohpm1 = 0.0;
	  if (p==1), ohpm1 = 1.0; end;
          for n=0:3,
              nm1 = n-1;
              % 0.0**n
              ohn = 0.0;
	      if (n==0), ohn = 1.0; end;
              % 0.0**(n-1)
              ohnm1 = 0.0,
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
                  ind = 16*p+4*n+m+1;
                  % m --> x
                  % n --> y
                  % p --> z
                  %-------------------------------------------------------------
                  %  VALUES
                  %-------------------------------------------------------------
                  A(1,ind) = ohm*ohn*ohp;
                  A(2,ind) = ohn*ohp;
                  A(3,ind) = ohm*ohp;
                  A(4,ind) = ohp;
                  A(5,ind) = ohm*ohn;
                  A(6,ind) = ohn;
                  A(7,ind) = ohm;
                  A(8,ind) = 1.0;
                  %-------------------------------------------------------------
                  %  Dx
                  %-------------------------------------------------------------
                  A(9,ind) = m*ohmm1*ohn*ohp;
                  A(10,ind) = m*ohn*ohp;
                  A(11,ind) = m*ohmm1*ohp;
                  A(12,ind) = m*ohp;
                  A(13,ind) = m*ohmm1*ohn;
                  A(14,ind) = m*ohn;
                  A(15,ind) = m*ohmm1;
                  A(16,ind) = m;
                  %-------------------------------------------------------------
                  %  Dy
                  %-------------------------------------------------------------
                  A(17,ind) = n*ohm*ohnm1*ohp;
                  A(18,ind) = n*ohnm1*ohp;
                  A(19,ind) = n*ohm*ohp;
                  A(20,ind) = n*ohp;
                  A(21,ind) = n*ohm*ohnm1;
                  A(22,ind) = n*ohnm1;
                  A(23,ind) = n*ohm;
                  A(24,ind) = n;
                  %-------------------------------------------------------------
                  %  Dz
                  %-------------------------------------------------------------
                  A(25,ind) = p*ohm*ohn*ohpm1;
                  A(26,ind) = p*ohn*ohpm1;
                  A(27,ind) = p*ohm*ohpm1;
                  A(28,ind) = p*ohpm1;
                  A(29,ind) = p*ohm*ohn;
                  A(30,ind) = p*ohn;
                  A(31,ind) = p*ohm;
                  A(32,ind) = p;
                  %-------------------------------------------------------------
                  %  DxDy
                  %-------------------------------------------------------------
                  f         = m*n;
                  A(33,ind) = f*ohmm1*ohnm1*ohp;
                  A(34,ind) = f*ohnm1*ohp;
                  A(35,ind) = f*ohmm1*ohp;
                  A(36,ind) = f*ohp;
                  A(37,ind) = f*ohmm1*ohnm1;
                  A(38,ind) = f*ohnm1;
                  A(39,ind) = f*ohmm1;
                  A(40,ind) = f;
                  %-------------------------------------------------------------
                  %  DxDz
                  %-------------------------------------------------------------
                  f         = m*p;
                  A(41,ind) = f*ohmm1*ohn*ohpm1;
                  A(42,ind) = f*ohn*ohpm1;
                  A(43,ind) = f*ohmm1*ohpm1;
                  A(44,ind) = f*ohpm1;
                  A(45,ind) = f*ohmm1*ohn;
                  A(46,ind) = f*ohn;
                  A(47,ind) = f*ohmm1;
                  A(48,ind) = f;
                  %-------------------------------------------------------------
                  %  DyDz
                  %-------------------------------------------------------------
                  f         = n*p;
                  A(49,ind) = f*ohm*ohnm1*ohpm1;
                  A(50,ind) = f*ohnm1*ohpm1;
                  A(51,ind) = f*ohm*ohpm1;
                  A(52,ind) = f*ohpm1;
                  A(53,ind) = f*ohm*ohnm1;
                  A(54,ind) = f*ohnm1;
                  A(55,ind) = f*ohm;
                  A(56,ind) = f;
                  %-------------------------------------------------------------
                  %  DxDyDz
                  %-------------------------------------------------------------
                  f         = m*n*p;
                  A(57,ind) = f*ohmm1*ohnm1*ohpm1;
                  A(58,ind) = f*ohnm1*ohpm1;
                  A(59,ind) = f*ohmm1*ohpm1;
                  A(60,ind) = f*ohpm1;
                  A(61,ind) = f*ohmm1*ohnm1;
                  A(62,ind) = f*ohnm1;
                  A(63,ind) = f*ohmm1;
                  A(64,ind) = f;
	      end;
	  end;
      end;

      [L,U,P] = lu(A);
      ALU = (L-eye(64))+U;
      ind = P*[1:1:64]';
      save 'gaga_ind.dat' ind -ASCII;
      ALUT = ALU';         % FORTRAN HAS COLUMN MAJOR
      save 'gaga_a.dat' ALUT -ASCII;

