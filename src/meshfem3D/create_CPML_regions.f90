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

  subroutine create_CPML_regions(nglob)

  use meshfem_par, only: ibool,prname, &
    nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML, &
    ADD_PML_AS_EXTRA_MESH_LAYERS, &
    SUPPRESS_UTM_PROJECTION,CREATE_VTK_FILES

  ! create the different regions of the mesh
  use constants, only: IMAIN,CUSTOM_REAL,SMALL_PERCENTAGE_TOLERANCE, &
    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ, &
    PI,TINYVAL,HUGEVAL,MAX_STRING_LEN,NDIM,myrank

  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use meshfem_par, only: nspec
  use create_meshfem_par, only: nodes_coords

  ! CPML
  use shared_parameters, only: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  implicit none

  integer,intent(inout):: nglob

  ! local parameters
  double precision :: xmin,xmax,ymin,ymax,zmin,zmax,limit
  double precision :: xmin_all,xmax_all,ymin_all,ymax_all,zmin_all,zmax_all
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
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1310')
  if (ier /= 0) stop 'Error allocating is_CPML array'
  ! initializes CPML elements
  is_CPML(:) = .false.
  nspec_CPML = 0

  ! checks if anything to do
  if (.not. PML_CONDITIONS) then
    ! dummy allocation
    allocate(CPML_to_spec(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1311')
    allocate(CPML_regions(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1312')
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

  ! PML as extra additional outer elements
  if (ADD_PML_AS_EXTRA_MESH_LAYERS) then
    ! create PML elements as additional outer layers
    call add_CPML_region_as_extra_layers(nglob)
    ! all done
    return
  endif

  ! here we assign PML elements within the existing mesh
  if (myrank == 0) then
    write(IMAIN,*) '  PML elements will be determined within mesh region and specified thicknesses'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! compute the min and max values of each coordinate
  xmin = minval(nodes_coords(:,1))
  xmax = maxval(nodes_coords(:,1))

  ymin = minval(nodes_coords(:,2))
  ymax = maxval(nodes_coords(:,2))

  zmin = minval(nodes_coords(:,3))
  zmax = maxval(nodes_coords(:,3))

  call min_all_all_dp(xmin,xmin_all)
  call min_all_all_dp(ymin,ymin_all)
  call min_all_all_dp(zmin,zmin_all)

  call max_all_all_dp(xmax,xmax_all)
  call max_all_all_dp(ymax,ymax_all)
  call max_all_all_dp(zmax,zmax_all)

  allocate(is_X_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1313')
  allocate(is_Y_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1314')
  allocate(is_Z_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1315')

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
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1316')
  allocate(CPML_regions(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1317')
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

!
!----------------------------------------------------------------------------------------------------
!

  subroutine add_CPML_region_as_extra_layers(nglob)

! PML as extra additional outer elements

  use meshfem_par, only: ibool,prname,CREATE_VTK_FILES,NEX_XI,NEX_ETA,NER, &
    UTM_X_MIN,UTM_X_MAX

  ! create the different regions of the mesh
  use constants, only: IMAIN,CUSTOM_REAL,SMALL_PERCENTAGE_TOLERANCE, &
    CPML_X_ONLY,CPML_Y_ONLY,CPML_Z_ONLY,CPML_XY_ONLY,CPML_XZ_ONLY,CPML_YZ_ONLY,CPML_XYZ, &
    PI,TINYVAL,HUGEVAL,MAX_STRING_LEN,NDIM,myrank, &
    GAUSSALPHA,GAUSSBETA

  use constants_meshfem, only: NGLLX_M,NGLLY_M,NGLLZ_M

  use meshfem_par, only: nspec
  use create_meshfem_par, only: nodes_coords,ispec_material_id, &
    iboun,iMPIcut_xi,iMPIcut_eta

  use shared_parameters, only: NGNOD

  ! CPML
  use shared_parameters, only: PML_CONDITIONS,PML_INSTEAD_OF_FREE_SURFACE

  use meshfem_par, only: nspec_CPML,is_CPML,CPML_to_spec,CPML_regions, &
    THICKNESS_OF_X_PML,THICKNESS_OF_Y_PML,THICKNESS_OF_Z_PML, &
    ADD_PML_AS_EXTRA_MESH_LAYERS,NUMBER_OF_PML_LAYERS_TO_ADD

  implicit none

  integer,intent(inout) :: nglob

  ! local parameters
  double precision :: xmin,xmax,ymin,ymax,zmin,zmax,limit
  double precision :: xmin_all,xmax_all,ymin_all,ymax_all,zmin_all,zmax_all
  double precision :: elem_size_x,elem_size_y,elem_size_z,elem_size

  logical, dimension(:), allocatable :: is_X_CPML,is_Y_CPML,is_Z_CPML
  logical, dimension(:), allocatable :: is_X_CPML_new,is_Y_CPML_new,is_Z_CPML_new

  integer :: nspec_CPML_total,nspec_total
  integer :: i1,i2,i3,i4,i5,i6,i7,i8
  integer :: ispec,ispec_CPML,ier
  character(len=MAX_STRING_LEN) :: filename

  ! for float comparisons
  !double precision, parameter :: TOL_EPS = 1.19d-7   ! machine precision for single precision (32-bit) floats
  !double precision :: rtol                           ! relative tolerance

  ! extra layers
  integer :: iloop_on_X_Y_Z_faces,iloop_on_min_face_then_max_face
  integer :: count_elem_faces_to_extend,count_elem_faces_to_extend_all
  integer :: npoin,nspec_all,nspec_new,npoin_new_max,npoin_new_real

  double precision :: xsize,ysize,zsize
  double precision :: value_min,value_max  !,value_size,sum_of_distances

  ! size of each PML element to add when they are added on the Xmin and Xmax faces, Ymin and Ymax faces, Zmin and/or Zmax faces
  double precision :: SIZE_OF_X_ELEMENT_TO_ADD,SIZE_OF_Y_ELEMENT_TO_ADD,SIZE_OF_Z_ELEMENT_TO_ADD

  logical :: ADD_ON_THE_XMIN_SURFACE,ADD_ON_THE_XMAX_SURFACE, &
             ADD_ON_THE_YMIN_SURFACE,ADD_ON_THE_YMAX_SURFACE, &
             ADD_ON_THE_ZMIN_SURFACE,ADD_ON_THE_ZMAX_SURFACE

  double precision, dimension(:), allocatable, target :: x,y,z
  double precision, dimension(:), pointer :: coord_to_use1,coord_to_use2,coord_to_use3

  double precision, dimension(:), allocatable :: x_new,y_new,z_new
  double precision, dimension(:,:), allocatable :: x_copy,y_copy,z_copy

  ! to check for negative Jacobians
  ! GLL points and weights of integration
  double precision :: xigll(NGLLX_M),yigll(NGLLY_M),zigll(NGLLZ_M), &
                      wxgll(NGLLX_M),wygll(NGLLY_M),wzgll(NGLLZ_M)
  ! 3D shape function derivatives,
  double precision :: shape3D(NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M), &
                      dershape3D(NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M)

  double precision, dimension(NGNOD) :: xelm,yelm,zelm
  double precision, dimension(NGNOD) :: x_tmp1,y_tmp1,z_tmp1,x_tmp2,y_tmp2,z_tmp2
  double precision :: jacobian

  integer, dimension(:,:,:,:), allocatable :: ibool_new
  integer, dimension(:), allocatable :: ispec_material_id_new

  ! boundary locator
  logical, dimension(:,:), allocatable :: iboun_new
  ! MPI cut-planes parameters along xi and along eta
  logical, dimension(:,:), allocatable :: iMPIcut_xi_new,iMPIcut_eta_new

  ! temporary array for local nodes
  integer, dimension(NGNOD) :: loc_node
  integer, dimension(NGNOD) :: anchor_iax,anchor_iay,anchor_iaz

  integer :: p1,p2,p3,p4,p5,p6,p7,p8,p9,ia,iglob
  integer :: factor_x,factor_y,factor_z
  integer :: iextend,elem_counter

  logical :: need_to_extend_this_element
  logical :: found_a_negative_Jacobian1,found_a_negative_Jacobian2

  ! tolerance for node detection on mesh boundaries
  double precision, parameter :: SMALL_RELATIVE_VALUE = 1.d-5

  ! for global indexing and mesh sorting
  integer :: ieoff,ilocnum,iglobnum
  ! variables for creating array ibool (some arrays also used for AVS or DX files)
  integer :: npointot
  integer, dimension(:), allocatable :: locval
  integer, dimension(:,:), allocatable :: ibool_tmp
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,zp

  ! (optional) updates absorbing boundary (iboun) array to move to outermost PML boundary elements
  ! if .false. then Stacey boundary stays at original mesh boundary
  logical, parameter :: MOVE_STACEY_ABSORBING_BOUNDARY = .false.
  ! sort and re-order new global nodes
  logical, parameter :: DO_MESH_SORTING = .true.

  ! safety check
  if (.not. PML_CONDITIONS) return
  if (.not. ADD_PML_AS_EXTRA_MESH_LAYERS) return

  if (myrank == 0) then
    write(IMAIN,*) '  adding PML layers as extra outer elements'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! check
  !if (ADD_PML_AS_EXTRA_MESH_LAYERS .and. NUMBER_OF_PML_LAYERS_TO_ADD == 0) &
  !  stop 'ADD_PML_AS_EXTRA_MESH_LAYERS requires NUMBER_OF_PML_LAYERS_TO_ADD > 0'

  ! initializes
  ADD_ON_THE_XMIN_SURFACE = .true.  ! left side
  ADD_ON_THE_XMAX_SURFACE = .true.  ! right
  ADD_ON_THE_YMIN_SURFACE = .true.  ! front
  ADD_ON_THE_YMAX_SURFACE = .true.  ! back
  ADD_ON_THE_ZMIN_SURFACE = .true.  ! bottom
  if (PML_INSTEAD_OF_FREE_SURFACE) then
    ADD_ON_THE_ZMAX_SURFACE = .true.  ! top
  else
    ADD_ON_THE_ZMAX_SURFACE = .false.
  endif

  ! sets up node addressing
  call hex_nodes_anchor_ijk_NGLL(NGNOD,anchor_iax,anchor_iay,anchor_iaz,NGLLX_M,NGLLY_M,NGLLZ_M)

  ! for jacobian testing
  ! set up coordinates of the Gauss-Lobatto-Legendre points
  call zwgljd(xigll,wxgll,NGLLX_M,GAUSSALPHA,GAUSSBETA)
  call zwgljd(yigll,wygll,NGLLY_M,GAUSSALPHA,GAUSSBETA)
  call zwgljd(zigll,wzgll,NGLLZ_M,GAUSSALPHA,GAUSSBETA)
  ! get the 3-D shape functions
  call get_shape3D(shape3D,dershape3D,xigll,yigll,zigll,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M)

  ! determine element sizes (all thicknesses given in m here)
  SIZE_OF_X_ELEMENT_TO_ADD = 0.d0
  SIZE_OF_Y_ELEMENT_TO_ADD = 0.d0
  SIZE_OF_Z_ELEMENT_TO_ADD = 0.d0
  if (NUMBER_OF_PML_LAYERS_TO_ADD > 0) then
    SIZE_OF_X_ELEMENT_TO_ADD = THICKNESS_OF_X_PML / NUMBER_OF_PML_LAYERS_TO_ADD
    SIZE_OF_Y_ELEMENT_TO_ADD = THICKNESS_OF_Y_PML / NUMBER_OF_PML_LAYERS_TO_ADD
    SIZE_OF_Z_ELEMENT_TO_ADD = THICKNESS_OF_Z_PML / NUMBER_OF_PML_LAYERS_TO_ADD
  endif

  if (myrank == 0) then
    write(IMAIN,*) '  number of PML layers to add = ',NUMBER_OF_PML_LAYERS_TO_ADD
    write(IMAIN,*) '  element type (NGNOD)        = ',NGNOD
    write(IMAIN,*)
    write(IMAIN,*) '  size of PML elements on X-faces to add = ',sngl(SIZE_OF_X_ELEMENT_TO_ADD),'m'
    write(IMAIN,*) '  size of PML elements on Y-faces to add = ',sngl(SIZE_OF_Y_ELEMENT_TO_ADD),'m'
    write(IMAIN,*) '  size of PML elements on Z-faces to add = ',sngl(SIZE_OF_Z_ELEMENT_TO_ADD),'m'
    write(IMAIN,*)
    if (ADD_ON_THE_ZMAX_SURFACE) then
      write(IMAIN,*) '  adding extra PML layers on top Z-face'
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! check
  if (NUMBER_OF_PML_LAYERS_TO_ADD == 0) then
    if (myrank == 0) then
      write(IMAIN,*) '***'
      write(IMAIN,*) '*** Warning: NUMBER_OF_PML_LAYERS_TO_ADD == 0, there are no PML layers to add (?)...'
      write(IMAIN,*) '***'
      write(IMAIN,*)
      write(IMAIN,*) 'no PML region'
      write(IMAIN,*)
      call flush_IMAIN()
    endif
    ! return here
    ! dummy allocation
    allocate(CPML_to_spec(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1311')
    allocate(CPML_regions(1),stat=ier)
    if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1312')
    if (ier /= 0) stop 'Error allocating dummy CPML arrays'
    ! all done, nothing to add
    return
  endif

  ! mesh point coordinates
  npoin = nglob

  allocate(x(npoin),y(npoin),z(npoin),stat=ier)
  if (ier /= 0) stop 'Error allocating x/y/z coordinates'
  x(:) = nodes_coords(:,1)
  y(:) = nodes_coords(:,2)
  z(:) = nodes_coords(:,3)

  ! initial total number of elements
  nspec_total = 0
  call sum_all_i(nspec,nspec_total)

  ! initial min and max values of each coordinate
  xmin = minval(x)
  xmax = maxval(x)

  ymin = minval(y)
  ymax = maxval(y)

  zmin = minval(z)
  zmax = maxval(z)

  ! determine global min/max
  call min_all_all_dp(xmin,xmin_all)
  call min_all_all_dp(ymin,ymin_all)
  call min_all_all_dp(zmin,zmin_all)
  xmin = xmin_all
  ymin = ymin_all
  zmin = zmin_all

  call max_all_all_dp(xmax,xmax_all)
  call max_all_all_dp(ymax,ymax_all)
  call max_all_all_dp(zmax,zmax_all)
  xmax = xmax_all
  ymax = ymax_all
  zmax = zmax_all

  ! mesh dimensions
  xsize = xmax - xmin
  ysize = ymax - ymin
  zsize = zmax - zmin

  ! estimanted element size
  elem_size_x = xsize / NEX_XI
  elem_size_y = ysize / NEX_ETA
  elem_size_z = ysize / NER

  ! minimum element size (for evaluation of point distances to mesh boundaries)
  elem_size = min(elem_size_x, elem_size_y)
  elem_size = min(elem_size, elem_size_z)

  if (myrank == 0) then
    write(IMAIN,*) '  initial mesh dimensions:'
    write(IMAIN,*) '    total number of mesh elements = ',nspec_total
    write(IMAIN,*)
    write(IMAIN,*) '    Xmin and Xmax of the mesh = ',sngl(xmin),sngl(xmax)
    write(IMAIN,*) '    Ymin and Ymax of the mesh = ',sngl(ymin),sngl(ymax)
    write(IMAIN,*) '    Zmin and Zmax of the mesh = ',sngl(zmin),sngl(zmax)
    write(IMAIN,*)
    write(IMAIN,*) '    size of the mesh along X = ',sngl(xsize),'m'
    write(IMAIN,*) '    size of the mesh along Y = ',sngl(ysize),'m'
    write(IMAIN,*) '    size of the mesh along Z = ',sngl(zsize),'m'
    write(IMAIN,*)
    write(IMAIN,*) '    estimated minimum mesh element size along X = ',sngl(elem_size_x),'m'
    write(IMAIN,*) '    estimated minimum mesh element size along Y = ',sngl(elem_size_y),'m'
    write(IMAIN,*) '    estimated minimum mesh element size along Z = ',sngl(elem_size_z),'m'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! flags for PML elements
  allocate(is_X_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1313')
  allocate(is_Y_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1314')
  allocate(is_Z_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating array 1315')
  is_X_CPML(:) = .false.
  is_Y_CPML(:) = .false.
  is_Z_CPML(:) = .false.

  ! C-PML element count
  nspec_CPML = 0

  ! loop on the three sets of faces to first add CPML elements along X, then along Y, then along Z
  do iloop_on_X_Y_Z_faces = 1,NDIM
    ! 1 is min face and 2 is max face (Xmin then Xmax, Ymin then Ymax, or Zmin then Zmax)
    do iloop_on_min_face_then_max_face = 1,2
      ! do not add a CPML layer on a given surface if the user asked not to
      if (iloop_on_X_Y_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_XMIN_SURFACE) cycle
      if (iloop_on_X_Y_Z_faces == 1 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_XMAX_SURFACE) cycle
      if (iloop_on_X_Y_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_YMIN_SURFACE) cycle
      if (iloop_on_X_Y_Z_faces == 2 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_YMAX_SURFACE) cycle
      if (iloop_on_X_Y_Z_faces == 3 .and. iloop_on_min_face_then_max_face == 1 .and. .not. ADD_ON_THE_ZMIN_SURFACE) cycle
      if (iloop_on_X_Y_Z_faces == 3 .and. iloop_on_min_face_then_max_face == 2 .and. .not. ADD_ON_THE_ZMAX_SURFACE) cycle

      if (myrank == 0) then
        write(IMAIN,*) '  ********************************************************************'
        if (iloop_on_X_Y_Z_faces == 1) then
          if (iloop_on_min_face_then_max_face == 1) then
            write(IMAIN,*) '  adding CPML elements along X-min faces of the existing mesh'
          else
            write(IMAIN,*) '  adding CPML elements along X-max faces of the existing mesh'
          endif
        else if (iloop_on_X_Y_Z_faces == 2) then
          if (iloop_on_min_face_then_max_face == 1) then
            write(IMAIN,*) '  adding CPML elements along Y-min faces of the existing mesh'
          else
            write(IMAIN,*) '  adding CPML elements along Y-max faces of the existing mesh'
          endif
        else if (iloop_on_X_Y_Z_faces == 3) then
          if (iloop_on_min_face_then_max_face == 1) then
            write(IMAIN,*) '  adding CPML elements along Z-min faces of the existing mesh'
          else
            write(IMAIN,*) '  adding CPML elements along Z-max faces of the existing mesh'
          endif
        else
          stop 'wrong index in loop on faces'
        endif
        write(IMAIN,*) '  ********************************************************************'
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      ! compute the current min and max values of each coordinate
      xmin = minval(x)
      xmax = maxval(x)

      ymin = minval(y)
      ymax = maxval(y)

      zmin = minval(z)
      zmax = maxval(z)

      ! determine global min/max
      call min_all_all_dp(xmin,xmin_all)
      call min_all_all_dp(ymin,ymin_all)
      call min_all_all_dp(zmin,zmin_all)
      xmin = xmin_all
      ymin = ymin_all
      zmin = zmin_all

      call max_all_all_dp(xmax,xmax_all)
      call max_all_all_dp(ymax,ymax_all)
      call max_all_all_dp(zmax,zmax_all)
      xmax = xmax_all
      ymax = ymax_all
      zmax = zmax_all

      ! mesh dimensions
      xsize = xmax - xmin
      ysize = ymax - ymin
      zsize = zmax - zmin

      if (myrank == 0) then
        write(IMAIN,*) '  current mesh dimensions:'
        write(IMAIN,*) '    Xmin and Xmax of the mesh = ',sngl(xmin),sngl(xmax)
        write(IMAIN,*) '    Ymin and Ymax of the mesh = ',sngl(ymin),sngl(ymax)
        write(IMAIN,*) '    Zmin and Zmax of the mesh = ',sngl(zmin),sngl(zmax)
        write(IMAIN,*)
        write(IMAIN,*) '    size of the mesh along X = ',sngl(xsize)
        write(IMAIN,*) '    size of the mesh along Y = ',sngl(ysize)
        write(IMAIN,*) '    size of the mesh along Z = ',sngl(zsize)
        write(IMAIN,*)
        call flush_IMAIN()
      endif

      ! determine coordinate system
      if (iloop_on_X_Y_Z_faces == 1) then
        ! Xmin or Xmax face
        value_min = xmin
        value_max = xmax
        coord_to_use1 => x  ! make coordinate array to use point to array x()
        coord_to_use2 => y
        coord_to_use3 => z
      else if (iloop_on_X_Y_Z_faces == 2) then
        ! Ymin or Ymax face
        value_min = ymin
        value_max = ymax
        coord_to_use1 => y  ! make coordinate array to use point to array y()
        coord_to_use2 => x
        coord_to_use3 => z
      else if (iloop_on_X_Y_Z_faces == 3) then
        ! Zmin or Zmax face
        value_min = zmin
        value_max = zmax
        coord_to_use1 => z  ! make coordinate array to use point to array z()
        coord_to_use2 => x
        coord_to_use3 => y
      else
        stop 'wrong index in loop on faces'
      endif

      ! check ibool
      !if (minval(ibool) /= 1) stop 'Error in minval(ibool) should be == 1'
      if (maxval(ibool) > nglob) stop 'Error in maxval(ibool) should be <= nglob'

      count_elem_faces_to_extend = 0
      !sum_of_distances = 0.d0

      ! loop on the whole mesh
      do ispec = 1,nspec
        ! we can use the 8 corners of the element only for the test here, even if the element is HEX27
        ! corner points
        i1 = ibool(1,1,1,ispec)
        i2 = ibool(NGLLX_M,1,1,ispec)
        i3 = ibool(NGLLX_M,NGLLY_M,1,ispec)
        i4 = ibool(1,NGLLY_M,1,ispec)
        i5 = ibool(1,1,NGLLZ_M,ispec)
        i6 = ibool(NGLLX_M,1,NGLLZ_M,ispec)
        i7 = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
        i8 = ibool(1,NGLLY_M,NGLLZ_M,ispec)

        if (iloop_on_min_face_then_max_face == 1) then
          ! min face
          ! detect elements belonging to the min face
          limit = value_min + elem_size * SMALL_RELATIVE_VALUE
          ! test face 1 (bottom)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
              coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
          endif
          ! test face 2 (top)
          if (coord_to_use1(i5) < limit .and. coord_to_use1(i6) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i5) - coord_to_use3(i6))**2 + (coord_to_use2(i5) - coord_to_use2(i6))**2)
          endif
          ! test face 3 (left)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit .and. &
              coord_to_use1(i8) < limit .and. coord_to_use1(i5) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i1) - coord_to_use3(i4))**2 + (coord_to_use2(i1) - coord_to_use2(i4))**2)
          endif
          ! test face 4 (right)
          if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i6) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i2) - coord_to_use3(i3))**2 + (coord_to_use2(i2) - coord_to_use2(i3))**2)
          endif
          ! test face 5 (front)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
              coord_to_use1(i6) < limit .and. coord_to_use1(i5) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
          endif
          ! test face 6 (back)
          if (coord_to_use1(i4) < limit .and. coord_to_use1(i3) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i4) - coord_to_use3(i3))**2 + (coord_to_use2(i4) - coord_to_use2(i3))**2)
          endif
        else
          ! max face
          ! detect elements belonging to the max face
          limit = value_max - elem_size * SMALL_RELATIVE_VALUE
          ! test face 1 (bottom)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
              coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !  sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
          endif
          ! test face 2 (top)
          if (coord_to_use1(i5) > limit .and. coord_to_use1(i6) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i5) - coord_to_use3(i6))**2 + (coord_to_use2(i5) - coord_to_use2(i6))**2)
          endif
          ! test face 3 (left)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit .and. &
              coord_to_use1(i8) > limit .and. coord_to_use1(i5) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i1) - coord_to_use3(i4))**2 + (coord_to_use2(i1) - coord_to_use2(i4))**2)
          endif
          ! test face 4 (right)
          if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i6) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i2) - coord_to_use3(i3))**2 + (coord_to_use2(i2) - coord_to_use2(i3))**2)
          endif
          ! test face 5 (front)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
              coord_to_use1(i6) > limit .and. coord_to_use1(i5) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i1) - coord_to_use3(i2))**2 + (coord_to_use2(i1) - coord_to_use2(i2))**2)
          endif
          ! test face 6 (back)
          if (coord_to_use1(i4) > limit .and. coord_to_use1(i3) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
            count_elem_faces_to_extend = count_elem_faces_to_extend + 1
            !sum_of_distances = sum_of_distances + &
            !    sqrt((coord_to_use3(i4) - coord_to_use3(i3))**2 + (coord_to_use2(i4) - coord_to_use2(i3))**2)
          endif
        endif
      enddo

      ! get total counts
      call sum_all_i(nspec,nspec_all)
      call sum_all_i(count_elem_faces_to_extend,count_elem_faces_to_extend_all)

      if (myrank == 0) then
        write(IMAIN,*) '  total number of elements in the mesh before extension = ',nspec_all
        write(IMAIN,*) '  number of element faces to extend  = ',count_elem_faces_to_extend_all
        write(IMAIN,*)
        write(IMAIN,*) '  adding a total number of PML elements to this face = ', &
                       count_elem_faces_to_extend_all * NUMBER_OF_PML_LAYERS_TO_ADD
        write(IMAIN,*)
        call flush_IMAIN()
        ! check total
        if (count_elem_faces_to_extend_all == 0) stop 'Error: number of element faces to extend detected is zero!'
      endif

      ! check if anything to do in this slice
      !if (count_elem_faces_to_extend == 0) cycle      ! make sure to have no more MPI collectives below this in this loop!!!

      ! we will add NUMBER_OF_LAYERS_TO_ADD to each of the element faces detected that need to be extended
      nspec_new = nspec + count_elem_faces_to_extend * NUMBER_OF_PML_LAYERS_TO_ADD

      !print *,'Total number of elements in the mesh after extension = ',nspec_new

      ! and each of these elements will have NGNOD points
      ! (some of them shared with other elements, but we will remove the multiples below, thus here it is a maximum
      npoin_new_max = npoin + count_elem_faces_to_extend * NGNOD * NUMBER_OF_PML_LAYERS_TO_ADD

      !mean_distance = sum_of_distances / dble(count_elem_faces_to_extend)
      !very_small_distance = mean_distance / 10000.d0
      !if (icompute_size == 1) print *,'Computed mean size of the elements to extend = ',mean_distance
      !print *

      ! check
      !if (minval(ibool) /= 1) stop 'Error in minval(ibool)'
      !if (maxval(ibool) > nglob) stop 'Error in maxval(ibool)'

      ! allocate a new set of elements, i.e. a new ibool()
      ! allocate arrays with the right size of the future extended mesh
      allocate(ispec_material_id_new(nspec_new), &
               ibool_new(NGLLX_M,NGLLY_M,NGLLZ_M,nspec_new),stat=ier)
      if (ier /= 0) stop 'Error allocating ibool_new array'
      ispec_material_id_new(:) = 0
      ibool_new(:,:,:,:) = 0
      ! copy the previous array values
      ispec_material_id_new(1:nspec) = ispec_material_id(1:nspec)
      ibool_new(:,:,:,1:nspec) = ibool(:,:,:,1:nspec)

      ! boundary & MPI
      allocate(iboun_new(6,nspec_new), &
               iMPIcut_xi_new(2,nspec_new), &
               iMPIcut_eta_new(2,nspec_new),stat=ier)
      if (ier /= 0) stop 'Error allocating iboun_new arrays'
      iboun_new(:,:) = .false.
      iMPIcut_xi_new(:,:) = .false.
      iMPIcut_eta_new(:,:) = .false.
      ! copy previous
      iboun_new(:,1:nspec) = iboun(:,:)
      iMPIcut_xi_new(:,1:nspec) = iMPIcut_xi(:,:)
      iMPIcut_eta_new(:,1:nspec) = iMPIcut_eta(:,:)

      ! C-PML arrays
      allocate(is_X_CPML_new(nspec_new), &
               is_Y_CPML_new(nspec_new), &
               is_Z_CPML_new(nspec_new), stat=ier)
      if (ier /= 0) stop 'Error allocating is_X_CPML_new array'
      is_X_CPML_new(:) = .false.; is_Y_CPML_new(:) = .false.; is_Z_CPML_new(:) = .false.
      ! copy previous
      is_X_CPML_new(1:nspec) = is_X_CPML(:)
      is_Y_CPML_new(1:nspec) = is_Y_CPML(:)
      is_Z_CPML_new(1:nspec) = is_Z_CPML(:)

      ! allocate a new set of points, with multiples
      allocate(x_new(npoin_new_max), &
               y_new(npoin_new_max), &
               z_new(npoin_new_max),stat=ier)
      if (ier /= 0) stop 'Error allocating x_new arrays'
      x_new(:) = 0.d0; y_new(:) = 0.d0; z_new(:) = 0.d0
      ! copy the original points into the new set
      x_new(1:npoin) = x(1:npoin)
      y_new(1:npoin) = y(1:npoin)
      z_new(1:npoin) = z(1:npoin)

      ! position after which to start to create the new elements
      elem_counter = nspec
      npoin_new_real = npoin

      ! loop on the whole original mesh
      do ispec = 1,nspec
        ! gets anchor nodes
        do ia = 1,NGNOD
          iglob = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
          loc_node(ia) = iglob
        enddo

  ! define topology of faces of cube for skewness
  ! (see hex_nodes.f90)
  !
  !           5 * -------- * 8
  !            /|         /|
  !         6 * -------- *7|
  !           |1* - - -  | * 4                 z
  !           |/         |/                   |
  !           * -------- *                    o--> y
  !          2           3                 x /

        ! corner points
        i1 = ibool(1,1,1,ispec)
        i2 = ibool(NGLLX_M,1,1,ispec)
        i3 = ibool(NGLLX_M,NGLLY_M,1,ispec)
        i4 = ibool(1,NGLLY_M,1,ispec)
        i5 = ibool(1,1,NGLLZ_M,ispec)
        i6 = ibool(NGLLX_M,1,NGLLZ_M,ispec)
        i7 = ibool(NGLLX_M,NGLLY_M,NGLLZ_M,ispec)
        i8 = ibool(1,NGLLY_M,NGLLZ_M,ispec)

        ! reset flag
        need_to_extend_this_element = .false.

        if (iloop_on_min_face_then_max_face == 1) then
          ! min face
          ! detect elements belonging to the min face
          limit = value_min + elem_size * SMALL_RELATIVE_VALUE
          ! test face 1 (bottom)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
              coord_to_use1(i3) < limit .and. coord_to_use1(i4) < limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i2
            p3 = i3
            p4 = i4
            if (NGNOD == 27) then
              p5 = loc_node(9)
              p6 = loc_node(10)
              p7 = loc_node(11)
              p8 = loc_node(12)
              p9 = loc_node(21)
            endif
          endif
          ! test face 2 (top)
          if (coord_to_use1(i5) < limit .and. coord_to_use1(i6) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
            need_to_extend_this_element = .true.
            p1 = i5
            p2 = i6
            p3 = i7
            p4 = i8
            if (NGNOD == 27) then
              p5 = loc_node(17)
              p6 = loc_node(18)
              p7 = loc_node(19)
              p8 = loc_node(20)
              p9 = loc_node(26)
            endif
          endif
          ! test face 3 (left)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i4) < limit .and. &
              coord_to_use1(i5) < limit .and. coord_to_use1(i8) < limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i4
            p3 = i8
            p4 = i5
            if (NGNOD == 27) then
              p5 = loc_node(12)
              p6 = loc_node(16)
              p7 = loc_node(20)
              p8 = loc_node(13)
              p9 = loc_node(25)
            endif
          endif
          ! test face 4 (right)
          if (coord_to_use1(i2) < limit .and. coord_to_use1(i3) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i6) < limit) then
            need_to_extend_this_element = .true.
            p1 = i2
            p2 = i3
            p3 = i7
            p4 = i6
            if (NGNOD == 27) then
              p5 = loc_node(10)
              p6 = loc_node(15)
              p7 = loc_node(18)
              p8 = loc_node(14)
              p9 = loc_node(23)
            endif
          endif
          ! test face 5 (front)
          if (coord_to_use1(i1) < limit .and. coord_to_use1(i2) < limit .and. &
              coord_to_use1(i6) < limit .and. coord_to_use1(i5) < limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i2
            p3 = i6
            p4 = i5
            if (NGNOD == 27) then
              p5 = loc_node(9)
              p6 = loc_node(14)
              p7 = loc_node(17)
              p8 = loc_node(13)
              p9 = loc_node(22)
            endif
          endif
          ! test face 6 (back)
          if (coord_to_use1(i4) < limit .and. coord_to_use1(i3) < limit .and. &
              coord_to_use1(i7) < limit .and. coord_to_use1(i8) < limit) then
            need_to_extend_this_element = .true.
            p1 = i4
            p2 = i3
            p3 = i7
            p4 = i8
            if (NGNOD == 27) then
              p5 = loc_node(11)
              p6 = loc_node(15)
              p7 = loc_node(19)
              p8 = loc_node(16)
              p9 = loc_node(24)
            endif
          endif
        else
          ! max face
          ! detect elements belonging to the max face
          limit = value_max - elem_size * SMALL_RELATIVE_VALUE
          ! test face 1 (bottom)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
              coord_to_use1(i3) > limit .and. coord_to_use1(i4) > limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i2
            p3 = i3
            p4 = i4
            if (NGNOD == 27) then
              p5 = loc_node(9)
              p6 = loc_node(10)
              p7 = loc_node(11)
              p8 = loc_node(12)
              p9 = loc_node(21)
            endif
          endif
          ! test face 2 (top)
          if (coord_to_use1(i5) > limit .and. coord_to_use1(i6) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
            need_to_extend_this_element = .true.
            p1 = i5
            p2 = i6
            p3 = i7
            p4 = i8
            if (NGNOD == 27) then
              p5 = loc_node(17)
              p6 = loc_node(18)
              p7 = loc_node(19)
              p8 = loc_node(20)
              p9 = loc_node(26)
            endif
          endif
          ! test face 3 (left)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i4) > limit .and. &
              coord_to_use1(i5) > limit .and. coord_to_use1(i8) > limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i4
            p3 = i8
            p4 = i5
            if (NGNOD == 27) then
              p5 = loc_node(12)
              p6 = loc_node(16)
              p7 = loc_node(20)
              p8 = loc_node(13)
              p9 = loc_node(25)
            endif
          endif
          ! test face 4 (right)
          if (coord_to_use1(i2) > limit .and. coord_to_use1(i3) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i6) > limit) then
            need_to_extend_this_element = .true.
            p1 = i2
            p2 = i3
            p3 = i7
            p4 = i6
            if (NGNOD == 27) then
              p5 = loc_node(10)
              p6 = loc_node(15)
              p7 = loc_node(18)
              p8 = loc_node(14)
              p9 = loc_node(23)
            endif
          endif
          ! test face 5 (front)
          if (coord_to_use1(i1) > limit .and. coord_to_use1(i2) > limit .and. &
              coord_to_use1(i6) > limit .and. coord_to_use1(i5) > limit) then
            need_to_extend_this_element = .true.
            p1 = i1
            p2 = i2
            p3 = i6
            p4 = i5
            if (NGNOD == 27) then
              p5 = loc_node(9)
              p6 = loc_node(14)
              p7 = loc_node(17)
              p8 = loc_node(13)
              p9 = loc_node(22)
            endif
          endif
          ! test face 6 (back)
          if (coord_to_use1(i4) > limit .and. coord_to_use1(i3) > limit .and. &
              coord_to_use1(i7) > limit .and. coord_to_use1(i8) > limit) then
            need_to_extend_this_element = .true.
            p1 = i4
            p2 = i3
            p3 = i7
            p4 = i8
            if (NGNOD == 27) then
              p5 = loc_node(11)
              p6 = loc_node(15)
              p7 = loc_node(19)
              p8 = loc_node(16)
              p9 = loc_node(24)
            endif
          endif
        endif

        if (need_to_extend_this_element) then
          ! create the NUMBER_OF_LAYERS_TO_ADD_IN_THIS_STEP new elements

          ! very important remark: it is OK to create duplicates of the mesh points in the loop below
          ! (i.e. not to tell the code that many of these points created are in fact shared between adjacent elements)
          ! because "xdecompose_mesh" will remove them automatically later on, thus no need to remove them here;
          ! this makes this PML mesh extrusion code much simpler to write.

          factor_x = 0
          factor_y = 0
          factor_z = 0

          if (iloop_on_X_Y_Z_faces == 1) then
            ! Xmin or Xmax
            if (iloop_on_min_face_then_max_face == 1) then
              ! min face
              factor_x = -1
            else
              ! max face
              factor_x = +1
            endif
          else if (iloop_on_X_Y_Z_faces == 2) then
            if (iloop_on_min_face_then_max_face == 1) then
              ! min face
              factor_y = -1
            else
              ! max face
              factor_y = +1
            endif
          else if (iloop_on_X_Y_Z_faces == 3) then
            if (iloop_on_min_face_then_max_face == 1) then
              ! min face
              factor_z = -1
            else
              ! max face
              factor_z = +1
            endif
          else
            stop 'wrong index in loop on faces'
          endif

          do iextend = 1,NUMBER_OF_PML_LAYERS_TO_ADD

            ! create a new element
            elem_counter = elem_counter + 1

            ! add to CPML
            nspec_CPML = nspec_CPML + 1

            ! assign CPML flags
            is_X_CPML_new(elem_counter) = is_X_CPML(ispec)
            is_Y_CPML_new(elem_counter) = is_Y_CPML(ispec)
            is_Z_CPML_new(elem_counter) = is_Z_CPML(ispec)
            if (iloop_on_X_Y_Z_faces == 1) then
              ! Xmin or Xmax
              is_X_CPML_new(elem_counter) = .true.
            else if (iloop_on_X_Y_Z_faces == 2) then
              ! Ymin or Ymax
              is_Y_CPML_new(elem_counter) = .true.
            else
              ! Zmin or Zmax
              is_Z_CPML_new(elem_counter) = .true.
            endif

            ! create the material property for the extended elements
            ispec_material_id_new(elem_counter) = ispec_material_id(ispec)

            ! boundary & MPI
            iboun_new(:,elem_counter) = .false. ! absorbing boundary for Stacey, can be left at original mesh boundary
            ! (optional) move absorbing boundary to outermost PML layer
            if (MOVE_STACEY_ABSORBING_BOUNDARY) then
              iboun(:,ispec) = .false.  ! reset as old boundary element is inside mesh now
              if (iextend == NUMBER_OF_PML_LAYERS_TO_ADD) then
                ! extended element is outer most element
                if (iloop_on_X_Y_Z_faces == 1) then
                  ! X-faces
                  if (iloop_on_min_face_then_max_face == 1) then
                    ! xmin
                    iboun_new(1,elem_counter) = .true.
                  else
                    ! xmax
                    iboun_new(2,elem_counter) = .true.
                  endif
                else if (iloop_on_X_Y_Z_faces == 2) then
                  ! Y-faces
                  if (iloop_on_min_face_then_max_face == 1) then
                    ! ymin
                    iboun_new(3,elem_counter) = .true.
                  else
                    ! ymax
                    iboun_new(4,elem_counter) = .true.
                  endif
                else
                  ! Z-faces
                  if (iloop_on_min_face_then_max_face == 1) then
                    ! zmin
                    iboun_new(5,elem_counter) = .true.
                  else
                    ! zmax
                    iboun_new(6,elem_counter) = .true.
                  endif
                endif
              else
                ! extended element is still inside mesh
                iboun_new(:,elem_counter) = .false.
              endif
            endif

            ! extend MPI-cut flags
            iMPIcut_xi_new(:,elem_counter) = iMPIcut_xi(:,ispec)
            iMPIcut_eta_new(:,elem_counter) = iMPIcut_eta(:,ispec)

            ! first create the normal element
            ! bottom face of HEX8
            x_tmp1(1) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp1(1) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp1(1) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp1(2) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp1(2) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp1(2) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp1(3) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp1(3) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp1(3) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp1(4) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp1(4) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp1(4) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            ! top face of HEX8
            x_tmp1(5) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp1(5) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp1(5) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp1(6) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp1(6) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp1(6) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp1(7) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp1(7) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp1(7) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp1(8) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp1(8) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp1(8) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            if (NGNOD == 27) then
              ! remaining points of bottom face of HEX27
              x_tmp1(9) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp1(9) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp1(9) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp1(10) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp1(10) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp1(10) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp1(11) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp1(11) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp1(11) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp1(12) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp1(12) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp1(12) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp1(21) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp1(21) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp1(21) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              ! remaining points of top face of HEX27
              x_tmp1(17) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp1(17) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp1(17) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp1(18) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp1(18) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp1(18) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp1(19) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp1(19) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp1(19) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp1(20) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp1(20) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp1(20) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp1(26) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp1(26) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp1(26) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              ! remaining points of middle cutplane (middle "face") of HEX27
              x_tmp1(13) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(13) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(13) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(14) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(14) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(14) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(15) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(15) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(15) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(16) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(16) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(16) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(22) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(22) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(22) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(23) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(23) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(23) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(24) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(24) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(24) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(25) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(25) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(25) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp1(27) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp1(27) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp1(27) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)
            endif

            ! then create the mirrored element
            ! bottom face of HEX8
            x_tmp2(1) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp2(1) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp2(1) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp2(2) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp2(2) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp2(2) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp2(3) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp2(3) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp2(3) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            x_tmp2(4) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
            y_tmp2(4) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
            z_tmp2(4) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

            ! top face of HEX8
            x_tmp2(5) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp2(5) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp2(5) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp2(6) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp2(6) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp2(6) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp2(7) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp2(7) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp2(7) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            x_tmp2(8) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
            y_tmp2(8) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
            z_tmp2(8) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

            if (NGNOD == 27) then
              ! remaining points of bottom face of HEX27
              x_tmp2(9) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp2(9) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp2(9) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp2(10) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp2(10) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp2(10) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp2(11) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp2(11) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp2(11) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp2(12) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp2(12) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp2(12) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              x_tmp2(21) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-1)
              y_tmp2(21) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-1)
              z_tmp2(21) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-1)

              ! remaining points of top face of HEX27
              x_tmp2(17) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp2(17) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp2(17) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp2(18) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp2(18) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp2(18) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp2(19) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp2(19) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp2(19) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp2(20) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp2(20) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp2(20) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              x_tmp2(26) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * iextend
              y_tmp2(26) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * iextend
              z_tmp2(26) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * iextend

              ! remaining points of middle cutplane (middle "face") of HEX27
              x_tmp2(13) = x(p1) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(13) = y(p1) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(13) = z(p1) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(14) = x(p2) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(14) = y(p2) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(14) = z(p2) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(15) = x(p3) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(15) = y(p3) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(15) = z(p3) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(16) = x(p4) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(16) = y(p4) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(16) = z(p4) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(22) = x(p5) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(22) = y(p5) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(22) = z(p5) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(23) = x(p6) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(23) = y(p6) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(23) = z(p6) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(24) = x(p7) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(24) = y(p7) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(24) = z(p7) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(25) = x(p8) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(25) = y(p8) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(25) = z(p8) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)

              x_tmp2(27) = x(p9) + factor_x * SIZE_OF_X_ELEMENT_TO_ADD * (iextend-0.5d0)
              y_tmp2(27) = y(p9) + factor_y * SIZE_OF_Y_ELEMENT_TO_ADD * (iextend-0.5d0)
              z_tmp2(27) = z(p9) + factor_z * SIZE_OF_Z_ELEMENT_TO_ADD * (iextend-0.5d0)
            endif

            ! now we need to test if the element created is flipped i.e. it has a negative Jacobian,
            ! and if so we will use the mirrored version of that element, which will then have a positive Jacobian

            ! check the element for a negative Jacobian
            do ia = 1,NGNOD
              xelm(ia) = x_tmp1(ia)
              yelm(ia) = y_tmp1(ia)
              zelm(ia) = z_tmp1(ia)
            enddo
            call check_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian1,NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M,jacobian)

            ! check the mirrored (i.e. flipped/swapped) element for a negative Jacobian
            ! either this one or the non-mirrored one above should be OK, and thus we will select it
            do ia = 1,NGNOD
              xelm(ia) = x_tmp2(ia)
              yelm(ia) = y_tmp2(ia)
              zelm(ia) = z_tmp2(ia)
            enddo
            call check_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian2,NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M,jacobian)

            ! this should never happen, it is just a safety test
            if (found_a_negative_jacobian1 .and. found_a_negative_jacobian2) then
              print *,'Error: rank ',myrank,'found a negative Jacobian that could not be fixed'
              print *,'       iloop_on_X_Y_Z_faces            :',iloop_on_X_Y_Z_faces
              print *,'       iloop_on_min_face_then_max_face :',iloop_on_min_face_then_max_face
              print *,'       iextend                         :',iextend
              print *,'       new element ',elem_counter
              print *,'       found negative jacobian: ',found_a_negative_jacobian1,found_a_negative_jacobian2
              print *,'       tmp1 : ',x_tmp1(1),y_tmp1(1),z_tmp1(1)
              print *,'       tmp2 : ',x_tmp2(1),y_tmp2(1),z_tmp2(1)
              stop 'Error: found a negative Jacobian that could not be fixed'
            endif

            ! this should also never happen, it is just a second safety test
            if (.not. found_a_negative_jacobian1 .and. .not. found_a_negative_jacobian2) &
              stop 'Error: both the element created and its mirrored version have a positive Jacobian!'

            ! if we have found that the original element has a negative Jacobian and its mirrored element is fine,
            ! swap the points so that we use that mirrored element in the final mesh saved to disk instead of the original one
            ! i.e. implement a mirror symmetry here
            if (found_a_negative_jacobian1) then
              do ia = 1,NGNOD
                npoin_new_real = npoin_new_real + 1
                !ibool_new(ia,elem_counter) = npoin_new_real
                ibool_new(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),elem_counter) = npoin_new_real
                x_new(npoin_new_real) = x_tmp2(ia)
                y_new(npoin_new_real) = y_tmp2(ia)
                z_new(npoin_new_real) = z_tmp2(ia)
              enddo
            else
              do ia = 1,NGNOD
                npoin_new_real = npoin_new_real + 1
                !ibool_new(ia,elem_counter) = npoin_new_real
                ibool_new(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),elem_counter) = npoin_new_real
                x_new(npoin_new_real) = x_tmp1(ia)
                y_new(npoin_new_real) = y_tmp1(ia)
                z_new(npoin_new_real) = z_tmp1(ia)
              enddo
            endif

            !debug
            !print *,'debug: iloop_on_X_Y_Z_faces',iloop_on_X_Y_Z_faces, &
            !        'iloop_on_min_face_then_max_face',iloop_on_min_face_then_max_face, &
            !        'iextend ',iextend,'added ',elem_counter,'nspec',nspec,nspec_new,'npoin_new_real',npoin_new_real, &
            !        'nglob',npoin,npoin_new_max, &
            !        'found jacobian',found_a_negative_jacobian1,found_a_negative_jacobian2
          enddo  ! iextend
        endif
      enddo

      ! checks
      !if (minval(ibool_new) /= 1) stop 'Error in minval(ibool_new)'
      if (maxval(ibool_new) > npoin_new_max) stop 'Error in maxval(ibool_new)'
      if (npoin_new_real /= npoin_new_max) stop 'Error in npoin_new_real'

      ! deallocate the original arrays
      deallocate(x,y,z)
      deallocate(ibool)
      deallocate(ispec_material_id)
      deallocate(iboun,iMPIcut_xi,iMPIcut_eta)
      deallocate(is_X_CPML,is_Y_CPML,is_Z_CPML)

      ! reallocate them with the new size
      allocate(x(npoin_new_real), &
               y(npoin_new_real), &
               z(npoin_new_real),stat=ier)
      if (ier /= 0) stop 'Error allocating new x,y,z,.. arrays'
      ! make the new ones become the old ones, to prepare for the next iteration of the two nested loops we are in,
      ! i.e. to make sure the next loop will extend the mesh from the new arrays rather than from the old ones
      x(:) = x_new(1:npoin_new_real)
      y(:) = y_new(1:npoin_new_real)
      z(:) = z_new(1:npoin_new_real)

      allocate(ispec_material_id(nspec_new), &
               ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec_new), &
               iboun(6,nspec_new), &
               iMPIcut_xi(2,nspec_new), &
               iMPIcut_eta(2,nspec_new), &
               is_X_CPML(nspec_new), &
               is_Y_CPML(nspec_new), &
               is_Z_CPML(nspec_new),stat=ier)
      if (ier /= 0) stop 'Error allocating new ispec_material_id,.. arrays'
      ispec_material_id(:) = ispec_material_id_new(:)
      ibool(:,:,:,:) = ibool_new(:,:,:,:)

      ! boundary & MPI
      iboun(:,:) = iboun_new(:,:)
      iMPIcut_xi(:,:) = iMPIcut_xi_new(:,:)
      iMPIcut_eta(:,:) = iMPIcut_eta_new(:,:)

      ! CPML flags
      is_X_CPML(:) = is_X_CPML_new(:)
      is_Y_CPML(:) = is_Y_CPML_new(:)
      is_Z_CPML(:) = is_Z_CPML_new(:)

      ! the new number of elements and points becomes the old one, for the same reason
      nspec = nspec_new
      npoin = npoin_new_real

      ! store as new nglob
      nglob = npoin

      ! check
      !if (minval(ibool) /= 1) stop 'Error in minval(ibool)'
      if (maxval(ibool) > nglob) stop 'Error in updated maxval(ibool)'

      ! deallocate the new ones, to make sure they can be allocated again in the next iteration of the nested loops we are in
      deallocate(x_new,y_new,z_new)
      deallocate(ibool_new)
      deallocate(ispec_material_id_new)
      deallocate(iboun_new,iMPIcut_xi_new,iMPIcut_eta_new)
      deallocate(is_X_CPML_new,is_Y_CPML_new,is_Z_CPML_new)

      ! new mesh dimensions
      xmin = minval(x)
      xmax = maxval(x)
      ymin = minval(y)
      ymax = maxval(y)
      zmin = minval(z)
      zmax = maxval(z)
      ! determine global min/max
      call min_all_all_dp(xmin,xmin_all)
      call min_all_all_dp(ymin,ymin_all)
      call min_all_all_dp(zmin,zmin_all)
      xmin = xmin_all
      ymin = ymin_all
      zmin = zmin_all
      call max_all_all_dp(xmax,xmax_all)
      call max_all_all_dp(ymax,ymax_all)
      call max_all_all_dp(zmax,zmax_all)
      xmax = xmax_all
      ymax = ymax_all
      zmax = zmax_all
      ! info
      if (myrank == 0) then
        write(IMAIN,*) '  new mesh dimensions:'
        write(IMAIN,*) '    Xmin and Xmax of the mesh = ',sngl(xmin),sngl(xmax)
        write(IMAIN,*) '    Ymin and Ymax of the mesh = ',sngl(ymin),sngl(ymax)
        write(IMAIN,*) '    Zmin and Zmax of the mesh = ',sngl(zmin),sngl(zmax)
        write(IMAIN,*)
        call flush_IMAIN()
      endif
    enddo  ! iloop_on_min_face_then_max_face
  enddo  ! iloop_on_X_Y_Z_faces

  ! we must remove all point multiples here using a fast sorting routine
  if (DO_MESH_SORTING) then
    if (myrank == 0) then
      write(IMAIN,*) '  beginning of multiple point removal based on sorting...'
      call flush_IMAIN()
    endif

    npointot = nspec * NGNOD

    allocate(ibool_tmp(NGNOD,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating ibool_tmp arrays'
    do ispec = 1,nspec
      do ia = 1,NGNOD
        ibool_tmp(ia,ispec) = ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec)
      enddo
    enddo

    allocate(locval(npointot), &
             ifseg(npointot), &
             xp(npointot), &
             yp(npointot), &
             zp(npointot),stat=ier)
    if (ier /= 0) stop 'Error allocating locval arrays'

    allocate(x_copy(NGNOD,nspec), &
             y_copy(NGNOD,nspec), &
             z_copy(NGNOD,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating x_copy arrays'

    ! puts x,y,z locations into 1D arrays
    do ispec = 1,nspec
      ieoff = NGNOD * (ispec-1)
      ilocnum = 0
      do ia = 1,NGNOD
        ilocnum = ilocnum + 1
        xp(ilocnum+ieoff) = x(ibool_tmp(ia,ispec))
        yp(ilocnum+ieoff) = y(ibool_tmp(ia,ispec))
        zp(ilocnum+ieoff) = z(ibool_tmp(ia,ispec))

        ! create a copy, since the sorting below will destroy the arrays
        x_copy(ia,ispec) = x(ibool_tmp(ia,ispec))
        y_copy(ia,ispec) = y(ibool_tmp(ia,ispec))
        z_copy(ia,ispec) = z(ibool_tmp(ia,ispec))
      enddo
    enddo

    ! gets ibool indexing from local (NGNOD points) to global points
    ! sorts xp,yp,zp in lexicographical order (increasing values)
    call get_global(npointot,xp,yp,zp,ibool_tmp,locval,ifseg,npoin,UTM_X_MIN,UTM_X_MAX)

    ! store as new nglob
    nglob = npoin

    if (myrank == 0) then
      write(IMAIN,*) '  done with multiple point removal based on sorting'
      write(IMAIN,*) '  found a total of ',npoin,' unique mesh points'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

    ! unique global point locations
    do ispec = 1, nspec
      do ia = 1,NGNOD
        iglobnum = ibool_tmp(ia,ispec)
        x(iglobnum) = x_copy(ia,ispec)
        y(iglobnum) = y_copy(ia,ispec)
        z(iglobnum) = z_copy(ia,ispec)
      enddo
    enddo

    ! re-allocate ibool
    deallocate(ibool)
    allocate(ibool(NGLLX_M,NGLLY_M,NGLLZ_M,nspec),stat=ier)
    if (ier /= 0) stop 'Error allocating sorted ibool'
    ! fill with sorted entries
    do ispec = 1,nspec
      do ia = 1,NGNOD
        ibool(anchor_iax(ia),anchor_iay(ia),anchor_iaz(ia),ispec) = ibool_tmp(ia,ispec)
      enddo
    enddo
  endif  ! DO_MESH_SORTING

  ! fill nodes_coords with added points
  deallocate(nodes_coords)

  ! note: nglob == npoin
  allocate(nodes_coords(nglob,NDIM),stat=ier)
  if (ier /= 0) stop 'Error allocating new nodes_coords array'
  nodes_coords(:,1) = x(1:nglob)
  nodes_coords(:,2) = y(1:nglob)
  nodes_coords(:,3) = z(1:nglob)

  ! free memory
  deallocate(x,y,z)

  ! final extended mesh min and max values of each coordinate
  xmin = minval(nodes_coords(:,1))
  xmax = maxval(nodes_coords(:,1))
  ymin = minval(nodes_coords(:,2))
  ymax = maxval(nodes_coords(:,2))
  zmin = minval(nodes_coords(:,3))
  zmax = maxval(nodes_coords(:,3))

  ! determine global min/max
  call min_all_all_dp(xmin,xmin_all)
  call min_all_all_dp(ymin,ymin_all)
  call min_all_all_dp(zmin,zmin_all)
  xmin = xmin_all
  ymin = ymin_all
  zmin = zmin_all
  call max_all_all_dp(xmax,xmax_all)
  call max_all_all_dp(ymax,ymax_all)
  call max_all_all_dp(zmax,zmax_all)
  xmax = xmax_all
  ymax = ymax_all
  zmax = zmax_all

  ! mesh dimensions
  xsize = xmax - xmin
  ysize = ymax - ymin
  zsize = zmax - zmin

  ! outputs total number of CPML elements
  nspec_total = 0
  call sum_all_i(nspec,nspec_total)

  nspec_CPML_total = 0
  call sum_all_i(nspec_CPML,nspec_CPML_total)

  call bcast_all_singlei(nspec_total)
  call bcast_all_singlei(nspec_CPML_total)

  if (myrank == 0) then
    write(IMAIN,*) '  resulting extended mesh dimensions:'
    write(IMAIN,*) '    total number of mesh elements = ',nspec_total
    write(IMAIN,*)
    write(IMAIN,*) '    Xmin and Xmax of the mesh = ',sngl(xmin),sngl(xmax)
    write(IMAIN,*) '    Ymin and Ymax of the mesh = ',sngl(ymin),sngl(ymax)
    write(IMAIN,*) '    Zmin and Zmax of the mesh = ',sngl(zmin),sngl(zmax)
    write(IMAIN,*)
    write(IMAIN,*) '    size of the mesh along X = ',sngl(xsize),'m'
    write(IMAIN,*) '    size of the mesh along Y = ',sngl(ysize),'m'
    write(IMAIN,*) '    size of the mesh along Z = ',sngl(zsize),'m'
    write(IMAIN,*)
    write(IMAIN,*) '  created a total of ',nspec_CPML_total,' unique CPML elements'
    write(IMAIN,*) '   (i.e., ',100. * nspec_CPML_total / real(nspec_total),'% of the mesh)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! re-new CPML array
  if (allocated(is_CPML)) deallocate(is_CPML)

  ! allocates C-PML arrays
  allocate(is_CPML(nspec),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating CPML_to_spec array')
  is_CPML(:) = .false.

  allocate(CPML_to_spec(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating CPML_to_spec array')
  CPML_to_spec(:) = 0

  allocate(CPML_regions(nspec_CPML),stat=ier)
  if (ier /= 0) call exit_MPI_without_rank('Error allocating CPML_regions array')
  if (ier /= 0) stop 'Error allocating CPML arrays'
  CPML_regions(:) = 0

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

  ! info
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) '  done adding C-PML elements as outer layers'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

contains

  subroutine check_jacobian(xelm,yelm,zelm,dershape3D,found_a_negative_jacobian,NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M,jacobian)

  implicit none

  integer,intent(in) :: NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M
  double precision, dimension(NGNOD),intent(in) :: xelm,yelm,zelm
  double precision,intent(in) :: dershape3D(NDIM,NGNOD,NGLLX_M,NGLLY_M,NGLLZ_M)
  logical,intent(out) :: found_a_negative_jacobian
  double precision,intent(out) :: jacobian
  ! local parameters
  integer :: i,j,k,ia
  double precision :: xxi,xeta,xgamma,yxi,yeta,ygamma,zxi,zeta,zgamma
  double precision, parameter :: ZERO = 0.d0

  found_a_negative_jacobian = .false.

  ! for this CPML mesh extrusion routine it is sufficient to test the 8 corners of each element to reduce the cost
  ! because we just want to detect if the element is flipped or not, and if so flip it back
  do k = 1,NGLLZ_M,NGLLZ_M-1
    do j = 1,NGLLY_M,NGLLY_M-1
      do i = 1,NGLLX_M,NGLLX_M-1

        xxi = ZERO
        xeta = ZERO
        xgamma = ZERO
        yxi = ZERO
        yeta = ZERO
        ygamma = ZERO
        zxi = ZERO
        zeta = ZERO
        zgamma = ZERO

        do ia = 1,NGNOD
          xxi = xxi + dershape3D(1,ia,i,j,k) * xelm(ia)
          xeta = xeta + dershape3D(2,ia,i,j,k) * xelm(ia)
          xgamma = xgamma + dershape3D(3,ia,i,j,k) * xelm(ia)

          yxi = yxi + dershape3D(1,ia,i,j,k) * yelm(ia)
          yeta = yeta + dershape3D(2,ia,i,j,k) * yelm(ia)
          ygamma = ygamma + dershape3D(3,ia,i,j,k) * yelm(ia)

          zxi = zxi + dershape3D(1,ia,i,j,k) * zelm(ia)
          zeta = zeta + dershape3D(2,ia,i,j,k) * zelm(ia)
          zgamma = zgamma + dershape3D(3,ia,i,j,k) * zelm(ia)
        enddo

        jacobian = xxi*(yeta*zgamma-ygamma*zeta) - xeta*(yxi*zgamma-ygamma*zxi) + xgamma*(yxi*zeta-yeta*zxi)

        ! check that the Jacobian transform is invertible, i.e. that the Jacobian never becomes negative or null
        if (jacobian <= ZERO) then
          found_a_negative_jacobian = .true.
          return
        endif

      enddo
    enddo
  enddo

  end subroutine check_jacobian

  end subroutine add_CPML_region_as_extra_layers
