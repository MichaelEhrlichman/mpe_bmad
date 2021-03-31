module noise_mod
  use bmad
  use bookkeeper_mod ! for intelligent bookkeeping
  
  implicit none

  integer, parameter :: spec_max_size = 4096
  
  type spec_data_struct
    real(rp) :: a   = -1.0
    real(rp) :: phi = 0.
    real(rp) :: f = 0.
    real(rp) :: df = 0.
  end type

  type nrm_data_struct
    real(rp) :: sigma = -1.0
  end type

  type ele_noise_struct_in
    character(200) :: ele_identifier = ''
    character(20) :: attrib_name = ''
    character(20) :: dist_type = 'fft'  !'fft', 'psd', 'nrm', 'nonuni_psd'
    character(20) :: error_type = 'add' ! must be add or mult or rel
    character(200) :: spec_data_file = 'none'
    character(10) :: phase_gang_element = 'none'
    type(spec_data_struct) :: spec_data(1:spec_max_size)
    integer :: n_turns_recalc = 0 ! how often to update magnets
    type(nrm_data_struct) :: nrm
    real(rp) :: dfreq = 0. ! step size between frequencies. a(1) is at frequency = start_freq
    real(rp) :: start_freq = 0. ! initial frequency
    logical :: random_phase_element = .false. 
    logical :: random_phase_component = .false. 
    !private
    integer n_loc
    character(20) :: bmad_custom_attribute = ''
    logical :: initialized = .false. 
    integer nspec !number of components in the spectra
  end type

  type, extends(ele_noise_struct_in) :: ele_noise_struct
    type (ele_pointer_struct), allocatable :: ma_eles(:)
  end type
    
contains
  ! subroutine init_ele_noise(ring, ele_noise)
  ! NOTE: this assumes ran_seed_put has already been called!
  subroutine init_ele_noise(ring, ele_noise)
    !use precision_def
    use random_mod
    
    implicit none 
    
    type(lat_struct), target :: ring
    type(ele_struct), pointer :: ele
    type(ele_noise_struct) :: ele_noise(:)
    type(all_pointer_struct) a_ptr
    integer :: ix, jx, kx,  n_ele_match = 0
    integer n_ele_noise
    real(rp) :: harvest
    logical err
    character(35) ele_set_str
    character(2) ix_str

    ! check if already initialized
    if (all(ele_noise(:)%initialized) .eqv. .true.) then
       write(*,*) 'init_ele_noise: ele_noise already initialized'
       return
    endif
      
    ! Count number of elements that this noise applies to.
    ! Store pointers to those elements in ma_eles
    n_ele_noise = 0
    do ix = 1, size(ele_noise)
      if (len(trim(ele_noise(ix)%attrib_name)) .eq. 0) exit
      n_ele_noise = n_ele_noise + 1
      call upcase_string(ele_noise(ix)%attrib_name)
      call downcase_string(ele_noise(ix)%error_type)
      call lat_ele_locator(ele_noise(ix)%ele_identifier, ring, ele_noise(ix)%ma_eles, ele_noise(ix)%n_loc) 
      if(ele_noise(ix)%n_loc .eq. 0) then
         write(*,'(a,a)') "Failed to find any elements matching identifier ", ele_noise(ix)%ele_identifier
      endif
    enddo

    ! Setup Bmad for ELE_PHASE_?? custom attribute
    ! This is a per-element phase offset for a given noise source
    ! For a given element, the same phase delta is applied to each component of the FFT spectrum
    do ix = 1, n_ele_noise
      write(ix_str,'(I2.2)') ix
      if( ele_noise(ix)%phase_gang_element == 'none' ) then
        ele_noise(ix)%bmad_custom_attribute = 'ELE_PHASE_'//ix_str
      else
        ele_noise(ix)%bmad_custom_attribute = 'ELE_PHASE_'//trim(ele_noise(ix)%phase_gang_element)
      endif
      ! the following call basically tests whether this attribute has been defined yet.  needed for ganged misalignments
      call set_ele_attribute(ele_noise(ix)%ma_eles(1)%ele, ele_noise(ix)%bmad_custom_attribute//'=0', err, .false.)
      if(err) then
        call set_custom_attribute_name(ele_noise(ix)%bmad_custom_attribute, err)
        call set_ele_attribute(ele_noise(ix)%ma_eles(1)%ele, ele_noise(ix)%bmad_custom_attribute//'=0', err, .true.)
        if(err) then
          write(*,'(4a)') "Failed to initialize custom attribute ", ele_noise(ix)%bmad_custom_attribute, &
                          " for element ", ring%ele(ix)%name
          stop
        endif
      endif
      do jx = 1, ele_noise(ix)%n_loc
        call set_ele_attribute(ele_noise(ix)%ma_eles(jx)%ele, ele_noise(ix)%bmad_custom_attribute//'=0', err, .false.)
      enddo
    enddo
    call lattice_bookkeeper(ring) ! probably unnecessary for ele_phase$, but good to be consistent. 

    if (all(len_trim(ele_noise(:)%attrib_name) .eq. 0)) then
      write(*,*) "All ele_noise(:)%attrib_name are zero. Skipping ele_noise_mod routines."
      return
    endif

    !count number of FFT or PSD components
    do jx = 1, n_ele_noise
      if( match_reg('fft|psd|nonuni_psd',ele_noise(jx)%dist_type)) then
        ele_noise(jx)%nspec = 0
        do ix=1, size(ele_noise(jx)%spec_data)
          if( ele_noise(jx)%spec_data(ix)%a .lt. 0.0 ) exit
          ele_noise(jx)%nspec = ele_noise(jx)%nspec + 1
        enddo
      endif
    enddo
        
    do jx = 1, n_ele_noise
      do ix = 1, ele_noise(jx)%n_loc
        ! For each noise source, for each element it is applied to
        ele => ele_noise(jx)%ma_eles(ix)%ele
        
        ! check if slave
        if (ele%slave_status .eq. super_slave$) then
          write(*,'(a,a)') "Is a slave: ", ele_noise(jx)%ele_identifier
          call err_exit()
        elseif (.not.(has_orientation_attributes(ele))) then
          write(*,'(a,a)') "Is not orientable: ", ele_noise(jx)%ele_identifier
          call err_exit()
        elseif (.not.attribute_free(ele, ele_noise(jx)%attrib_name, err_print_flag = .false.)) then
          write(*,'(4a)') "Attribute not free. ele identifier: ", ele_noise(jx)%ele_identifier, &
                             "attribute name: ", ele_noise(jx)%attrib_name
          call err_exit()
        endif

        ! if random_phase_element for the source is true, populate ELE_PHASE_?? with the random phase
        if(ele_noise(jx)%random_phase_element) then
          call ran_uniform(harvest) ! generates flat distribution on [0,1]
          write(ele_set_str,'(a,es14.4)') ele_noise(jx)%bmad_custom_attribute//'=', harvest*twopi
          call set_ele_attribute(ele, ele_set_str, err, .true.)
          if(err) then
            write(*,*) "bad"
            stop
          endif
        endif
      enddo
    enddo

    do jx = 1, n_ele_noise
      if( match_reg('nonuni_psd',ele_noise(jx)%dist_type)) then
        do ix=2,ele_noise(jx)%nspec-1
          ele_noise(jx)%spec_data(ix)%df = ( ele_noise(jx)%spec_data(ix+1)%f - ele_noise(jx)%spec_data(ix-1)%f ) / 2.0d0
        enddo
      endif
    enddo

    ! revisit: each component of a PSD is given a random phase
    do jx = 1, n_ele_noise
      if( match_reg('fft|psd|nonuni_psd',ele_noise(jx)%dist_type)) then
        if( ele_noise(jx)%random_phase_component ) then
          do ix = 1, ele_noise(jx)%nspec
            call ran_uniform(harvest) ! generates flat distribution on [0,1]
            ele_noise(jx)%spec_data(ix)%phi = twopi * harvest
          enddo
        endif
      endif
    enddo
    call lattice_bookkeeper(ring) ! probably unnecessary for ele_phase$, but good to be consistent. 

    ele_noise(:)%initialized = .true.

    !do jx=1,n_ele_noise
    !  do ix=1,ele_noise(jx)%n_loc
    !    ele => ele_noise(jx)%ma_eles(ix)%ele
    !    write(*,*) "FOO custom phi: ", value_of_attribute(ele,ele_noise(jx)%bmad_custom_attribute,err) 
    !  enddo
    !  do ix=1,n_ele_noise
    !    write(*,*) "FOO phi: ", ele_noise(jx)%spec_data(ix)%phi
    !  enddo
    !enddo
  end subroutine

  !---------------------------------------------------------
  ! subroutine apply_ele_noise(ring, ring_ideal, ele_noise, turn_number)
  ! subroutine to apply periodic noise to a lattice
  !---------------------------------------------------------
  subroutine apply_ele_noise(ring, ring_ideal, ele_noise, turn_number)
    !use precision_def

    implicit none 

    type(lat_struct), target :: ring
    type(lat_struct) :: ring_ideal
    type(ele_struct), pointer :: ele, ele_ideal
    type(ele_noise_struct) :: ele_noise(:)
    type (all_pointer_struct) :: a_ptr, a_ptr_ideal
    integer turn_number
    integer :: ix, jx, kx
    real(rp) :: net_kick 
    real(rp) :: circ_time = 0., dfreq = 0.
    real(rp) ideal_val
    logical err
    character(100) ele_set_str

    circ_time = ring%ele(ring%n_ele_track)%s / c_light

    do jx = 1, size(ele_noise)
      if (mod(turn_number,ele_noise(jx)%n_turns_recalc) .ne. 0) cycle

      ! now cycle through elements to apply the updated kicks
      do ix = 1, ele_noise(jx)%n_loc
        ele => ele_noise(jx)%ma_eles(ix)%ele
        ele_ideal => ring_ideal%ele(ele%ix_ele)

        if (ele%slave_status .eq. super_slave$) cycle
        if (.not. (has_orientation_attributes(ele) .or. &
             attribute_free(ele, ele_noise(jx)%attrib_name, err_print_flag = .false.))) cycle

        call calc_noise_kick(ele_noise(jx), ele, turn_number, circ_time, net_kick)

        if(match_reg('G_ERR',ele_noise(jx)%attrib_name)) then
          ideal_val = value_of_attribute(ele_ideal, 'G' ,err_print_flag = .true.)
        else
          ideal_val = value_of_attribute(ele_ideal, ele_noise(jx)%attrib_name,err_print_flag = .true.)
        endif

        if( match_reg('fft|psd|nonuni_psd',ele_noise(jx)%dist_type)) then
          if (match_reg('add',ele_noise(jx)%error_type)) then
            write(ele_set_str,'(a,a,es16.9)') ele_noise(jx)%attrib_name,'=', ideal_val + net_kick
          elseif (match_reg('mult',ele_noise(jx)%error_type)) then
            write(ele_set_str,'(a,a,es16.9)') ele_noise(jx)%attrib_name,'=', ideal_val * (1.0d0 + net_kick)
          elseif (match_reg('rel',ele_noise(jx)%error_type)) then
            write(ele_set_str,'(a,a,es16.9)') ele_noise(jx)%attrib_name,'=', ideal_val * net_kick
          else
             write(*,*) "Error applying noise.  Check attrib_name. Stopping here..."
             call err_exit()
          endif
        elseif(match_reg('nrm',ele_noise(jx)%dist_type)) then
            write(ele_set_str,'(a,a,es14.4)') ele_noise(jx)%attrib_name,'=', ideal_val * (1.0d0 + net_kick)
        else
          write(*,*) "Unable to match dist_type: ", ele_noise(jx)%dist_type
          call err_exit()
        endif

        call set_ele_attribute(ele, ele_set_str, err, .true.)
        if(err) write(*,*) "Error"
      enddo
    enddo 
    call lattice_bookkeeper(ring)
  end subroutine
  
  !###################################################
  subroutine calc_noise_kick(ele_noise, ele, turn_number, circ_time, net_kick)
    use precision_def
    
    implicit none 
    
    type(ele_noise_struct) ele_noise
    type(ele_struct) :: ele
    real(rp) net_kick, circ_time, dfreq
    real(rp) harvest
    real(rp) per_ele_phase
    real(rp) freq(ele_noise%nspec)
    logical err
    integer i, turn_number

    dfreq = ele_noise%dfreq

    if(ele_noise%random_phase_element) then
      per_ele_phase = value_of_attribute(ele,ele_noise%bmad_custom_attribute,err)
      if(err) then
        write(*,*) "That was bad."
        stop
      endif
    else
      per_ele_phase = 0.0d0
    endif

    if(match_reg('fft|psd',ele_noise%dist_type)) then
      freq(:) = ele_noise%start_freq + (/(i, i=0,ele_noise%nspec-1)/)*dfreq
      if(match_reg('fft',ele_noise%dist_type)) then
        net_kick = sum(ele_noise%spec_data(1:ele_noise%nspec)%a* &
                       sin(twopi*freq(:)*(circ_time*turn_number) + ele_noise%spec_data(1:ele_noise%nspec)%phi + per_ele_phase))
      elseif(match_reg('psd',ele_noise%dist_type)) then
        net_kick = sum(sqrt(2.0d0*ele_noise%spec_data(1:ele_noise%nspec)%a*dfreq)* &
                       sin(twopi*freq(:)*(circ_time*turn_number) + ele_noise%spec_data(1:ele_noise%nspec)%phi + per_ele_phase))
      endif
    elseif(match_reg('nonuni_psd',ele_noise%dist_type)) then
      associate(sd => ele_noise%spec_data(1:ele_noise%nspec))
        net_kick = sum(sqrt(2.0d0*sd%a*sd%df) * sin(twopi*sd%f*(circ_time*turn_number) + sd%phi + per_ele_phase))
      end associate
    elseif(match_reg('nrm',ele_noise%dist_type)) then
      call ran_gauss(harvest)
      net_kick = harvest*ele_noise%nrm%sigma
    endif

  end subroutine
end module
