!=====================================================================
!
!                          S p e c f e m 3 D
!                          -----------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

  subroutine create_CPML_regions(nspec,nglob,nodes_coords)

  use meshfem_par, only: ibool,prname, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML, &
    SUPPRESS_UTM_PROJECTION,CREATE_VTK_FILES

  ! create the different regions of the mesh
  use constants, only: IMAIN,CUSTOM_REAL,SMALL_PERCENTAGE_TOLERANCE, &
    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ, &
    PI,TINYVAL,HUGEVAL,MAX_STRING_LEN,NDIM,myrank

  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  ! CPML
  use shared_parameters, only: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  implicit none

  integer,intent(in):: nspec,nglob
  double precision, dimension(nglob,NDIM), intent(in) :: nodes_coords

  ! local parameters
  real(kind=CUSTOM_REAL) :: xmin,xmax,ymin,ymax,zmin,zmax,limit
  real(kind=CUSTOM_REAL) :: xmin_all,xmax_all,ymin_all,ymax_all,zmin_all,zmax_all
  double precision :: elem_size_x_min,elem_size_x_max
  double precision :: elem_size_y_min,elem_size_y_max
  double precision :: elem_size_z_min,elem_size_z_max
  double precision :: elem_size,elem_size_all
  logical, dimension(:), allocatable :: is_X_CPML,is_Y_CPML,is_Z_CPML

  integer :: nspec_CPML_total,nspec_total
  integer :: i1,i2,i3,i4,i5,i6,i7,i8
  integer :: ispec,ispec_CPML
  integer :: ier
  character(len=MAX_STRING_LEN) :: filename

  ! conversion factor from degrees to m (1 degree = 6371.d0 * PI/180 = 111.1949 km) and vice versa
  double precision, parameter :: DEGREES_TO_METERS = 6371000.d0 * PI/180.d0
  double precision, parameter :: METERS_TO_DEGREES = 1.d0 / (6371000.d0 * PI/180.d0)

  ! for float comparisons
  double precision, parameter :: TOL_EPS = 1.19d-7   ! machine precision for single precision (32-bit) floats
  double precision :: rtol                           ! relative tolerance

  ! CPML allocation
  allocate(is_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1310')
  if (ier /= 0) stop 'Error allocating is_CPML array'
  ! initializes CPML elements
  is_CPML(:) = .false.
  nspec_CPML = 0

  ! checks if anything to do
  if (.not. PML_CONDITIONS) then
    ! dummy allocation
    allocate(CPML_to_spec(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1311')
    allocate(CPML_regions(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('error allocating array 1312')
    if (ier /= 0) stop 'Error allocating dummy CPML arrays'

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) 'no PML region'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! nothing to do anymore
    return
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'creating PML region'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (SUPPRESS_UTM_PROJECTION) then
    ! input lat/lon given as X/Y directly in m
    if (myrank == 0) then
      write(IMAIN,*) '  THICKNESS_OF_X_PML (in m) = ',sngl(THICKNESS_OF_X_PML)
      write(IMAIN,*) '  THICKNESS_OF_Y_PML (in m) = ',sngl(THICKNESS_OF_Y_PML)
      write(IMAIN,*) '  THICKNESS_OF_Z_PML (in m) = ',sngl(THICKNESS_OF_Z_PML)
      write(IMAIN,*)
    endif
  else
    ! using UTM projection, all locations and widths given in input file are in degree
    if (myrank == 0) then
      write(IMAIN,*) '  THICKNESS_OF_X_PML (in degree) = ',sngl(THICKNESS_OF_X_PML)
      write(IMAIN,*) '  THICKNESS_OF_Y_PML (in degree) = ',sngl(THICKNESS_OF_Y_PML)
      write(IMAIN,*) '  THICKNESS_OF_Z_PML (in degree) = ',sngl(THICKNESS_OF_Z_PML)
      write(IMAIN,*)
    endif
    ! checks
    if (THICKNESS_OF_X_PML >= 360.d0) stop 'Error invalid degree value for THICKNESS_OF_X_PML'
    if (THICKNESS_OF_Y_PML >= 360.d0) stop 'Error invalid degree value for THICKNESS_OF_Y_PML'
    if (THICKNESS_OF_Z_PML >= 360.d0) stop 'Error invalid degree value for THICKNESS_OF_Z_PML'

    ! converts thickness to m
    THICKNESS_OF_X_PML = THICKNESS_OF_X_PML * DEGREES_TO_METERS
    THICKNESS_OF_Y_PML = THICKNESS_OF_Y_PML * DEGREES_TO_METERS
    THICKNESS_OF_Z_PML = THICKNESS_OF_Z_PML * DEGREES_TO_METERS
    if (myrank == 0) then
      write(IMAIN,*) '  using UTM projection, thickness converted to meters:'
      write(IMAIN,*) '  THICKNESS_OF_X_PML (in m) = ',sngl(THICKNESS_OF_X_PML)
      write(IMAIN,*) '  THICKNESS_OF_Y_PML (in m) = ',sngl(THICKNESS_OF_Y_PML)
      write(IMAIN,*) '  THICKNESS_OF_Z_PML (in m) = ',sngl(THICKNESS_OF_Z_PML)
      write(IMAIN,*)
    endif
  endif

  ! compute the min and max values of each coordinate
  xmin = minval(nodes_coords(:,1))
  xmax = maxval(nodes_coords(:,1))

  ymin = minval(nodes_coords(:,2))
  ymax = maxval(nodes_coords(:,2))

  zmin = minval(nodes_coords(:,3))
  zmax = maxval(nodes_coords(:,3))

  call min_all_all_cr(xmin,xmin_all)
  call min_all_all_cr(ymin,ymin_all)
  call min_all_all_cr(zmin,zmin_all)

  call max_all_all_cr(xmax,xmax_all)
  call max_all_all_cr(ymax,ymax_all)
  call max_all_all_cr(zmax,zmax_all)

  allocate(is_X_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1313')
  allocate(is_Y_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1314')
  allocate(is_Z_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1315')

  is_X_CPML(:) = .false.
  is_Y_CPML(:) = .false.
  is_Z_CPML(:) = .false.

  ! stats
  elem_size_x_min = HUGEVAL
  elem_size_x_max = 0.d0
  elem_size_y_min = HUGEVAL
  elem_size_y_max = 0.d0
  elem_size_z_min = HUGEVAL
  elem_size_z_max = 0.d0

  do ispec = 1,nspec
    ! corner points
    i1 = ibool(1,1,1,ispec)
    i2 = ibool(NGLLX_M,1,1,ispec)
    i3 = ibool(NGLLX_M,NGLLY_M,1,ispec)
    i4 = ibool(1,NGLLY_M,1,ispec)
    i5 = ibool(1,1,NGLLZ_M,ispec)
    i6 = ibool(NGLLX_M,1,NGLLZ_M,ispec)
    i7 = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
    i8 = ibool(1,NGLLY_M,NGLLZ_M,ispec)

    ! note: if user THICKNESS_OF_X_PML/THICKNESS_OF_Y_PML/THICKNESS_OF_Z_PML was specified as > 0,
    !       then we add at least the outer most elements of the mesh, even if these elements are larger
    !       than the specified thickness.
    !       this avoids having no PML elements at the mesh boundary when doubling layers are used
    !       and element sizes vary with depth.

    ! mesh sides along X
    if (THICKNESS_OF_X_PML > 0.d0) then
      ! Xmin CPML
      limit = xmin_all + THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
      if ( nodes_coords(i1,1) < limit .and. nodes_coords(i2,1) < limit .and. &
           nodes_coords(i3,1) < limit .and. nodes_coords(i4,1) < limit .and. &
           nodes_coords(i5,1) < limit .and. nodes_coords(i6,1) < limit .and. &
           nodes_coords(i7,1) < limit .and. nodes_coords(i8,1) < limit) then
        is_X_CPML(ispec) = .true.
      endif
      ! Xmax CPML
      limit = xmax_all - THICKNESS_OF_X_PML * SMALL_PERCENTAGE_TOLERANCE
      if ( nodes_coords(i1,1) > limit .and. nodes_coords(i2,1) > limit .and. &
           nodes_coords(i3,1) > limit .and. nodes_coords(i4,1) > limit .and. &
           nodes_coords(i5,1) > limit .and. nodes_coords(i6,1) > limit .and. &
           nodes_coords(i7,1) > limit .and. nodes_coords(i8,1) > limit) then
        is_X_CPML(ispec) = .true.
      endif

      ! check if on X-min side (i==1)
      if (abs(xmin_all) > 0.0_CUSTOM_REAL) then
        rtol = abs(xmin_all) * TOL_EPS  ! float comparison tolerance
      else
        rtol = TOL_EPS                  ! machine precision
      endif
      if ( abs(nodes_coords(i1,1) - xmin_all) < rtol &
          .or. abs(nodes_coords(i4,1) - xmin_all) < rtol &
          .or. abs(nodes_coords(i5,1) - xmin_all) < rtol &
          .or. abs(nodes_coords(i8,1) - xmin_all) < rtol ) then
        ! element has corner(s) on X-min side
        is_X_CPML(ispec) = .true.
      endif
      ! check if on X-max side (i==NGLLX)
      if (abs(xmax_all) > 0.0_CUSTOM_REAL) then
        rtol = abs(xmax_all) * TOL_EPS  ! float comparison tolerance
      else
        rtol = TOL_EPS                  ! machine precision
      endif
      if ( abs(nodes_coords(i2,1) - xmax_all) < rtol &
          .or. abs(nodes_coords(i3,1) - xmax_all) < rtol &
          .or. abs(nodes_coords(i6,1) - xmax_all) < rtol &
          .or. abs(nodes_coords(i7,1) - xmax_all) < rtol ) then
        ! element has corner(s) on X-max side
        is_X_CPML(ispec) = .true.
      endif
    endif

    ! mesh sides along Y
    if (THICKNESS_OF_Y_PML > 0.d0) then
      ! Ymin CPML
      limit = ymin_all + THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
      if ( nodes_coords(i1,2) < limit .and. nodes_coords(i2,2) < limit .and. &
           nodes_coords(i3,2) < limit .and. nodes_coords(i4,2) < limit .and. &
           nodes_coords(i5,2) < limit .and. nodes_coords(i6,2) < limit .and. &
           nodes_coords(i7,2) < limit .and. nodes_coords(i8,2) < limit) then
        is_Y_CPML(ispec) = .true.
      endif
      ! Ymax CPML
      limit = ymax_all - THICKNESS_OF_Y_PML * SMALL_PERCENTAGE_TOLERANCE
      if ( nodes_coords(i1,2) > limit .and. nodes_coords(i2,2) > limit .and. &
           nodes_coords(i3,2) > limit .and. nodes_coords(i4,2) > limit .and. &
           nodes_coords(i5,2) > limit .and. nodes_coords(i6,2) > limit .and. &
           nodes_coords(i7,2) > limit .and. nodes_coords(i8,2) > limit) then
        is_Y_CPML(ispec) = .true.
      endif

      ! check if on Y-min side (j==1)
      if (abs(ymin_all) > 0.0_CUSTOM_REAL) then
        rtol = abs(ymin_all) * TOL_EPS  ! float comparison tolerance
      else
        rtol = TOL_EPS                  ! machine precision
      endif
      if ( abs(nodes_coords(i1,2) - ymin_all) < rtol &
          .or. abs(nodes_coords(i2,2) - ymin_all) < rtol &
          .or. abs(nodes_coords(i5,2) - ymin_all) < rtol &
          .or. abs(nodes_coords(i6,2) - ymin_all) < rtol ) then
        ! element has corner(s) on Y-min side
        is_Y_CPML(ispec) = .true.
      endif
      ! check if on Y-max side (j==NGLLY)
      if (abs(ymax_all) > 0.0_CUSTOM_REAL) then
        rtol = abs(ymax_all) * TOL_EPS  ! float comparison tolerance
      else
        rtol = TOL_EPS                  ! machine precision
      endif
      if ( abs(nodes_coords(i3,2) - ymax_all) < rtol &
          .or. abs(nodes_coords(i4,2) - ymax_all) < rtol &
          .or. abs(nodes_coords(i7,2) - ymax_all) < rtol &
          .or. abs(nodes_coords(i8,2) - ymax_all) < rtol ) then
        ! element has corner(s) on Y-max side
        is_Y_CPML(ispec) = .true.
      endif
    endif

    ! mesh sides along Z
    if (THICKNESS_OF_Z_PML > 0.d0) then
      ! Zmin CPML
      limit = zmin_all + THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
      if ( nodes_coords(i1,3) < limit .and. nodes_coords(i2,3) < limit .and. &
           nodes_coords(i3,3) < limit .and. nodes_coords(i4,3) < limit .and. &
           nodes_coords(i5,3) < limit .and. nodes_coords(i6,3) < limit .and. &
           nodes_coords(i7,3) < limit .and. nodes_coords(i8,3) < limit) then
         is_Z_CPML(ispec) = .true.
      endif
      ! check if on Z-min side (k==1)
      if (abs(zmin_all) > 0.0_CUSTOM_REAL) then
        rtol = abs(zmin_all) * TOL_EPS  ! float comparison tolerance
      else
        rtol = TOL_EPS                  ! machine precision
      endif
      if ( abs(nodes_coords(i1,3) - zmin_all) < rtol &
          .or. abs(nodes_coords(i2,3) - zmin_all) < rtol &
          .or. abs(nodes_coords(i3,3) - zmin_all) < rtol &
          .or. abs(nodes_coords(i4,3) - zmin_all) < rtol ) then
        ! element has corner(s) on X-min side
        is_Z_CPML(ispec) = .true.
      endif

      if (PML_INSTEAD_OF_FREE_SURFACE) then
        ! Zmax CPML
        limit = zmax_all - THICKNESS_OF_Z_PML * SMALL_PERCENTAGE_TOLERANCE
        if (  nodes_coords(i1,3) > limit .and. nodes_coords(i2,3) > limit .and. &
              nodes_coords(i3,3) > limit .and. nodes_coords(i4,3) > limit .and. &
              nodes_coords(i5,3) > limit .and. nodes_coords(i6,3) > limit .and. &
              nodes_coords(i7,3) > limit .and. nodes_coords(i8,3) > limit) then
          is_Z_CPML(ispec) = .true.
        endif
        ! check if on Z-max side (k==NGLLZ)
        if (abs(zmax_all) > 0.0_CUSTOM_REAL) then
          rtol = abs(zmax_all) * TOL_EPS  ! float comparison tolerance
        else
          rtol = TOL_EPS                  ! machine precision
        endif
        if ( abs(nodes_coords(i5,3) - zmax_all) < rtol &
            .or. abs(nodes_coords(i6,3) - zmax_all) < rtol &
            .or. abs(nodes_coords(i7,3) - zmax_all) < rtol &
            .or. abs(nodes_coords(i8,3) - zmax_all) < rtol ) then
          ! element has corner(s) on Z-max side
          is_Z_CPML(ispec) = .true.
        endif
      endif
    endif

    ! total counter
    if (is_X_CPML(ispec) .or. is_Y_CPML(ispec) .or. is_Z_CPML(ispec)) nspec_CPML = nspec_CPML + 1

    ! stats
    if (is_X_CPML(ispec)) then
      ! element size in X-direction
      elem_size = abs(nodes_coords(i1,1)-nodes_coords(i2,1))                     ! (1,1,1)-(NGLLX,1,1)
      elem_size = max(elem_size,abs(nodes_coords(i4,1)-nodes_coords(i3,1)))      ! (1,NGLLY,1)-(NGLLX,NGLLY,1)
      elem_size = max(elem_size,abs(nodes_coords(i5,1)-nodes_coords(i6,1)))      ! (1,1,NGLLZ)-(NGLLX,1,NGLLZ)
      elem_size = max(elem_size,abs(nodes_coords(i8,1)-nodes_coords(i7,1)))      ! (1,NGLLY,NGLLZ)-(NGLLX,NGLLY,NGLLZ)
      ! overall min/max
      elem_size_x_min = min(elem_size_x_min,elem_size)
      elem_size_x_max = max(elem_size_x_max,elem_size)
    endif
    if (is_Y_CPML(ispec)) then
      ! element size in y-direction
      elem_size = abs(nodes_coords(i1,2)-nodes_coords(i4,2))                     ! (1,1,1)-(1,NGLLY,1)
      elem_size = max(elem_size,abs(nodes_coords(i2,2)-nodes_coords(i3,2)))      ! (NGLLX,1,1)-(NGLLX,NGLLY,1)
      elem_size = max(elem_size,abs(nodes_coords(i5,2)-nodes_coords(i8,2)))      ! (1,1,NGLLZ)-(1,NGLLY,NGLLZ)
      elem_size = max(elem_size,abs(nodes_coords(i6,2)-nodes_coords(i7,2)))      ! (NGLLX,1,NGLLZ)-(NGLLX,NGLLY,NGLLZ)
      ! overall min/max
      elem_size_y_min = min(elem_size_y_min,elem_size)
      elem_size_y_max = max(elem_size_y_max,elem_size)
    endif
    if (is_Z_CPML(ispec)) then
      ! element size in y-direction
      elem_size = abs(nodes_coords(i1,3)-nodes_coords(i5,3))                     ! (1,1,1)-(1,1,NGLLZ)
      elem_size = max(elem_size,abs(nodes_coords(i2,3)-nodes_coords(i6,3)))      ! (NGLLX,1,1)-(NGLLX,1,NGLLZ)
      elem_size = max(elem_size,abs(nodes_coords(i4,3)-nodes_coords(i8,3)))      ! (1,NGLLY,1)-(1,NGLLY,NGLLZ)
      elem_size = max(elem_size,abs(nodes_coords(i3,3)-nodes_coords(i7,3)))      ! (NGLLX,NGLLY,1)-(NGLLX,NGLLY,NGLLZ)
      ! overall min/max
      elem_size_z_min = min(elem_size_z_min,elem_size)
      elem_size_z_max = max(elem_size_z_max,elem_size)
    endif
  enddo

  ! outputs total number of CPML elements
  nspec_total = 0
  call sum_all_i(nspec,nspec_total)

  nspec_CPML_total = 0
  call sum_all_i(nspec_CPML,nspec_CPML_total)

  call bcast_all_singlei(nspec_total)
  call bcast_all_singlei(nspec_CPML_total)

  if (myrank == 0) then
    write(IMAIN,*) '  Created a total of ',nspec_CPML_total,' unique CPML elements'
    write(IMAIN,*) '   (i.e., ',100. * nspec_CPML_total / real(nspec_total),'% of the mesh)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! outputs stats of CPML elements
  ! X side elements
  call min_all_all_dp(elem_size_x_min,elem_size_all)
  elem_size_x_min = elem_size_all
  call max_all_all_dp(elem_size_x_max,elem_size_all)
  elem_size_x_max = elem_size_all
  ! Y side elements
  call min_all_all_dp(elem_size_y_min,elem_size_all)
  elem_size_y_min = elem_size_all
  call max_all_all_dp(elem_size_y_max,elem_size_all)
  elem_size_y_max = elem_size_all
  ! Z side elements
  call min_all_all_dp(elem_size_z_min,elem_size_all)
  elem_size_z_min = elem_size_all
  call max_all_all_dp(elem_size_z_max,elem_size_all)
  elem_size_z_max = elem_size_all

  if (myrank == 0) then
    write(IMAIN,*) '  found PML elements: min/max element size along X sides (in m)      = ', &
                   sngl(elem_size_x_min),'/',sngl(elem_size_x_max)
    write(IMAIN,*) '                      min/max element size along Y sides (in m)      = ', &
                   sngl(elem_size_y_min),'/',sngl(elem_size_y_max)
    write(IMAIN,*) '                      min/max element size along Z sides (in m)      = ', &
                   sngl(elem_size_z_min),'/',sngl(elem_size_z_max)
    write(IMAIN,*)
    ! converts to degrees
    if (.not. SUPPRESS_UTM_PROJECTION) then
      write(IMAIN,'(a,f8.5,a,f8.5)') '                       min/max element size along X sides (in degree) = ',&
                   sngl(elem_size_x_min*METERS_TO_DEGREES),'/',sngl(elem_size_x_max*METERS_TO_DEGREES)
      write(IMAIN,'(a,f8.5,a,f8.5)') '                       min/max element size along Y sides (in degree) = ',&
                   sngl(elem_size_y_min*METERS_TO_DEGREES),'/',sngl(elem_size_y_max*METERS_TO_DEGREES)
      write(IMAIN,'(a,f8.5,a,f8.5)') '                       min/max element size along Z sides (in degree) = ',&
                   sngl(elem_size_z_min*METERS_TO_DEGREES),'/',sngl(elem_size_z_max*METERS_TO_DEGREES)
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! allocates arrays
  allocate(CPML_to_spec(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1316')
  allocate(CPML_regions(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('error allocating array 1317')
  if (ier /= 0) stop 'Error allocating CPML arrays'

  ispec_CPML = 0
  do ispec = 1,nspec
    if (is_X_CPML(ispec) .and. is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_XYZ
       is_CPML(ispec) = .true.
    else if (is_Y_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_YZ_ONLY
       is_CPML(ispec) = .true.

    else if (is_X_CPML(ispec) .and. is_Z_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_XZ_ONLY
       is_CPML(ispec) = .true.

    else if (is_X_CPML(ispec) .and. is_Y_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_XY_ONLY
       is_CPML(ispec) = .true.

    else if (is_Z_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_Z_ONLY
       is_CPML(ispec) = .true.

    else if (is_Y_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_Y_ONLY
       is_CPML(ispec) = .true.

    else if (is_X_CPML(ispec)) then
       ispec_CPML = ispec_CPML+1
       CPML_to_spec(ispec_CPML) = ispec
       CPML_regions(ispec_CPML) = CPML_X_ONLY
       is_CPML(ispec) = .true.
   endif
  enddo

  ! checks
  if (ispec_CPML /= nspec_CPML) stop 'Error number of CPML element is not consistent'

  ! file output
  if (CREATE_VTK_FILES) then
    ! vtk file output
    filename = prname(1:len_trim(prname))//'is_CPML.vtk'
    if (myrank == 0) then
      write(IMAIN,*) '  saving VTK file: ',trim(filename)
    endif
    call write_VTK_data_elem_i_meshfemCPML(nglob,nspec,NGLLX_M,nodes_coords,ibool, &
                                           nspec_CPML,CPML_regions,is_CPML,filename)
  endif

  end subroutine create_CPML_regions
