    DO k=kmin,kmax 
      DO j=jmin,jmax 
        DO i=imin,imax 
          !---------------------------------------------------------------------
          ! First check if the point is inside or outside the solid
          !---------------------------------------------------------------------
          vecT(1) = REAL(i-1+(ptcset(ilevel,ipatch)%sistr(1,isubl)-1),mk)*dx+ptcset(1,1)%min(1)
          vecT(2) = REAL(j-1+(ptcset(ilevel,ipatch)%sistr(2,isubl)-1),mk)*dy+ptcset(1,1)%min(2)
          vecT(3) = REAL(k-1+(ptcset(ilevel,ipatch)%sistr(3,isubl)-1),mk)*dz+ptcset(1,1)%min(3)

          intersections = 0
          inout         = 0
          DO tr=1,stlt(1)%tri_count
            !check if triangle is aligned with the point in the bounds_direction
            !perhaps add elseif verbose a more thorough test
            ! maybe change to some finite value
            !implement warning if stl points are too close to grid coordinates
            IF (stlt(1)%tri_norm(tr,__z).NE.0.0_mk) THEN 
              vecX(1) = vecT(__x) - stlt(1)%tri_base(tr,__x)
              vecX(2) = vecT(__y) - stlt(1)%tri_base(tr,__y)
              vecX(3) = vecT(__z) - stlt(1)%tri_base(tr,__z)

              !beware of vecX index
              udotx = stlt(1)%tri_vecu(tr,__x)*vecX(1) + &
                    & stlt(1)%tri_vecu(tr,__y)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecu(tr,__z)*vecX(3)
              vdotx = stlt(1)%tri_vecv(tr,__x)*vecX(1) + &
                    & stlt(1)%tri_vecv(tr,__y)*vecX(2) !@+ &
                    !@& stlt(1)%tri_vecv(tr,__z)*vecX(3)

              !get barycentric coordinates
              a = (udotx*stlt(1)%tri_vdotv2d(tr)-vdotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              b = (vdotx*stlt(1)%tri_udotu2d(tr)-udotx*stlt(1)%tri_udotv2d(tr))*stlt(1)%tri_denom2d(tr)
              c = (1.0_mk-a-b)
              IF (a.GE.0.0_mk .AND. b.GE.0.0_mk.AND.c.GE.0.0_mk) THEN
                dist = a*(stlt(1)%tri_base(tr,__z)+stlt(1)%tri_vecu(tr,__z)) + &
                     & b*(stlt(1)%tri_base(tr,__z)+stlt(1)%tri_vecv(tr,__z)) + &
                     & c* stlt(1)%tri_base(tr,__z)
                IF (dist .LT. vecT(__z)) THEN
                  intersections = intersections + 1
                  IF (stlopt_check_intersections) THEN
                    IF (stlt(1)%tri_norm(tr,__z).GT.0.0_mk) THEN
                      inout = inout + 1
                    ELSE
                      inout = inout - 1
                    ENDIF
                  ENDIF
                ENDIF
              ENDIF
            ENDIF
          ENDDO !tr

          minimumdist = diag
          DO tr=1,stlt(1)%tri_count
            vecS=stlt(1)%tri_base(tr,:)-vecT
            !----------------------------------------------------------------
            !          x3          x0 is the projection of T onto the
            !          /\             triangle plane
            ! x0      /  \         x1 is used as triangle base point
            !        /    \        vecU = x2 - x1
            !       /      \       vecV = x3 - x1
            !      /________\      vecW = x3 - x2
            !     x1        x2
            !
            !    x0   x1                 T is the point to be evaluated
            !    - - - ======== - - -    === is the triangle (in its plane)
            !     |   /                  vecP = x0 - x1
            !vecP |  / vecS              vecS = x1 - T
            !     | /
            !     |/
            !     T
            !
            ! Meassure distance to triangle.
            ! This is not the vecP from figure 6.2 in thesis
            ! a and b are being reused. Initial definition is in thesis.
            ! tri_norm   contains normal vector of triangle plane
            !----------------------------------------------------------------
            vecP = stlt(1)%tri_norm(tr,:) * &
            & ( stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
            & stlt(1)%tri_norm(tr,3)*vecS(3)) - vecS

            Pdotu = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecu(tr,3)
            Pdotv = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
            & vecP(3)*stlt(1)%tri_vecv(tr,3)

            a = (Pdotu*stlt(1)%tri_vdotv(tr) - Pdotv*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)
            b = (Pdotv*stlt(1)%tri_udotu(tr) - Pdotu*stlt(1)%tri_udotv(tr)) * stlt(1)%tri_denom(tr)

            IF ((a .GE. 0.0_mk) .AND. (b .GE. 0.0_mk) .AND. & 
            & ((a+b) .LE. 1.0_mk)) THEN
              !get distance normal to triangle
              a = (stlt(1)%tri_norm(tr,1)*vecS(1) + stlt(1)%tri_norm(tr,2)*vecS(2) + &
              & stlt(1)%tri_norm(tr,3)*vecS(3))**2
            ELSE
              !get distance to edge or corner
              a = vecP(1)*stlt(1)%tri_vecu(tr,1) + vecP(2)*stlt(1)%tri_vecu(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecu(tr,3)
              b = vecP(1)*stlt(1)%tri_vecv(tr,1) + vecP(2)*stlt(1)%tri_vecv(tr,2) + &
                & vecP(3)*stlt(1)%tri_vecv(tr,3)
              c = (vecP(1)-stlt(1)%tri_vecu(tr,1))*stlt(1)%tri_vecw(tr,1) + &
                & (vecP(2)-stlt(1)%tri_vecu(tr,2))*stlt(1)%tri_vecw(tr,2) + &
                & (vecP(3)-stlt(1)%tri_vecu(tr,3))*stlt(1)%tri_vecw(tr,3)

              IF (a .GT. stlt(1)%tri_udotu(tr)) THEN !point is closest to x2
                !a = tri_udotu(tr)
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                a = SUM(vece**2)
              ELSEIF (a .LT. 0.0_mk) THEN !point is closest to x1
                !a = 0
                !vece=-vecS
                a = SUM(vecS**2)
              ELSE !point is closest to vecU
                vece=-vecS-stlt(1)%tri_vecu(tr,:)*a/stlt(1)%tri_udotu(tr)
                a = SUM(vece**2)
                !this may not be the cheapest way to do this
              ENDIF 

              IF (b .GT. stlt(1)%tri_vdotv(tr)) THEN !point is closest to x3
                !b = tri_vdotv(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                b = SUM(vece**2)
              ELSEIF (b .LT. 0.0_mk) THEN !point is closest to x1
                !b = 0
                !vece=-vecS
                b = SUM(vecS**2)
              ELSE !point is closest to vecV
                vece=-vecS-stlt(1)%tri_vecv(tr,:)*b/stlt(1)%tri_vdotv(tr)
                b = SUM(vece**2)
              ENDIF 

              IF (c .GT. stlt(1)%tri_wdotw(tr)) THEN !point is closest to x3
                !c = tri_wdotw(tr)
                vece=-vecS-stlt(1)%tri_vecv(tr,:)
                c = SUM(vece**2)
              ELSEIF (c .LT. 0.0_mk) THEN !point is closest to x2
                !c = 0
                vece=-vecS-stlt(1)%tri_vecu(tr,:)
                c = SUM(vecS**2)
              ELSE !point is closest to vecW
                vece=-vecS-stlt(1)%tri_vecu(tr,:)-stlt(1)%tri_vecw(tr,:)*c/stlt(1)%tri_wdotw(tr)
                c = SUM(vece**2)
              ENDIF

              a = min(a,b)
              a = min(a,c)

            ENDIF

            minimumdist = min(minimumdist,a)
          END DO !tr triangles

          !!chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)
          minimumdist = sqrt(minimumdist)!JTR cleanup
          IF (stlopt_check_intersections) THEN
            IF (inout.LT.0) THEN
              minimumdist = sign(minimumdist,-1.0_mk)
            ENDIF
          ELSE IF(MOD(intersections,2) .NE. 0) THEN
            minimumdist = sign(minimumdist,-1.0_mk)
          ENDIF

          !The multiplication with epsilon can be done to the step function limits instead
          !Could be done either in this routine or when the limits are determined (better). 
          !That would save a multiplication of epsilon, but would in the long run require 
          !multiple epsilon. Leave it for now
          chif(ilevel,ipatch)%fld(i,j,k,isub) = naga_stepfunc1(minimumdist*epsilon)

        END DO !inner array index
      END DO !middle array index
    END DO !outer array index
