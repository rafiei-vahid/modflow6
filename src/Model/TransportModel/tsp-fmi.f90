module TspFmiModule

  use KindModule, only: DP, I4B, LGP
  use ConstantsModule, only: DONE, DZERO, DHALF, LINELENGTH, LENBUDTXT, &
                             LENPACKAGENAME, LENVARNAME
  use SimModule, only: store_error, count_errors, store_error_filename
  use SimVariablesModule, only: errmsg
  use FlowModelInterfaceModule, only: FlowModelInterfaceType
  use BaseDisModule, only: DisBaseType
  use ListModule, only: ListType
  use BudgetFileReaderModule, only: BudgetFileReaderType
  use HeadFileReaderModule, only: HeadFileReaderType
  use PackageBudgetModule, only: PackageBudgetType
  use BudgetObjectModule, only: BudgetObjectType, budgetobject_cr_bfr
  use MatrixBaseModule

  implicit none
  private
  public :: TspFmiType
  public :: fmi_cr

  character(len=LENPACKAGENAME) :: text = '    GWTFMI'

  integer(I4B), parameter :: NBDITEMS = 2
  character(len=LENBUDTXT), dimension(NBDITEMS) :: budtxt
  data budtxt/'      FLOW-ERROR', ' FLOW-CORRECTION'/

  type :: DataAdvancedPackageType
    real(DP), dimension(:), contiguous, pointer :: concpack => null()
    real(DP), dimension(:), contiguous, pointer :: qmfrommvr => null()
  end type

  type :: BudObjPtrArray
    type(BudgetObjectType), pointer :: ptr
  end type BudObjPtrArray

  !> E v3 (rafiei-vahid fork): the final destinations of one E-inactive cell's inflowing mass
  type :: RouteDestType
    integer(I4B), dimension(:), allocatable :: t !< 1 = active receiver cell, 2 = sink at an E-inactive cell
    integer(I4B), dimension(:), allocatable :: n !< the receiver node (1) or the E-inactive node holding the sink (2)
    real(DP), dimension(:), allocatable :: f !< fraction of the cell's inflowing mass reaching it
    integer(I4B), dimension(:), allocatable :: p !< receivers: flowja position of the face (receiver, adjacent E-inactive cell)
  end type RouteDestType

  type, extends(FlowModelInterfaceType) :: TspFmiType

    integer(I4B), dimension(:), pointer, contiguous :: iatp => null() !< advanced transport package applied to gwfpackages
    integer(I4B), pointer :: iflowerr => null() !< add the flow error correction
    real(DP), dimension(:), pointer, contiguous :: flowcorrect => null() !< mass flow correction
    real(DP), pointer :: eqnsclfac => null() !< governing equation scale factor; =1. for solute; =rhow*cpw for energy
    type(DataAdvancedPackageType), &
      dimension(:), pointer, contiguous :: datp => null()
    type(BudObjPtrArray), dimension(:), allocatable :: aptbudobj !< flow budget objects for the advanced packages
    integer(I4B) :: izerosatinact = 0 !< patch E (rafiei-vahid fork): 1 = a cell at zero flow saturation is inactive for transport
    integer(I4B), dimension(:), allocatable :: idrynew !< patch E: 1 where set_active_status deactivated the cell THIS step
    ! -- E v3 (rafiei-vahid fork): advected mass follows the water through E-inactive cells (DRY_CELL_ROUTE_MASS)
    integer(I4B) :: iroute = 0 !< 1 = route the mass (requires izerosatinact)
    integer(I4B), dimension(:), allocatable :: rstat !< per node: 0 not E-inactive, 1 routable E-inactive, 2 isolated E-inactive
    integer(I4B), dimension(:), allocatable :: rpos !< per node: position in rinode (0 if not E-inactive)
    integer(I4B) :: nri = 0 !< number of E-inactive cells this step
    integer(I4B), dimension(:), allocatable :: rinode !< the E-inactive cells
    integer(I4B), dimension(:), allocatable :: rdia !< per rinode position: start of its final-destination list (nri + 1)
    integer(I4B), dimension(:), allocatable :: rdtype !< 1 = active receiver, 2 = sink
    integer(I4B), dimension(:), allocatable :: rdnode !< receiver node (1) or the E-inactive node holding the sink (2)
    real(DP), dimension(:), allocatable :: rdfrac !< fraction of the inflowing mass reaching the destination
    integer(I4B), dimension(:), allocatable :: rdpos !< receivers: flowja position in the receiver's row of its face with
                                                     !! the adjacent E-inactive cell (the routed mass enters there)
    real(DP), dimension(:), allocatable :: raptlag !< per node: advanced-package mass into an E-inactive cell, fc iterate
    integer(I4B) :: naptc = 0 !< converged advanced-package connections into E-inactive cells (this step)
    integer(I4B), dimension(:), allocatable :: raptf, raptm !< feature and cell of each
    real(DP), dimension(:), allocatable :: raptw, raptq !< water rate and mass rate (Q x C_feature, converged) of each
    integer(I4B) :: nro = 0 !< ROUTED list length
    integer(I4B) :: nrs = 0 !< ROUTED-SINK list length
    integer(I4B) :: nrw = 0 !< WITHHELD-ISOL list length
    integer(I4B), dimension(:), allocatable :: ronode, roid2, rsnode, rsid2, rwnode, rwid2
    real(DP), dimension(:), allocatable :: roq, rsq, rwq
    real(DP) :: rsumd_fc = DZERO !< sum of the deliveries the last outer iteration applied
    real(DP) :: rcncin = DZERO !< routed mass leaving constant-concentration senders (the boundary supplies it)
    real(DP) :: rcncout = DZERO !< routed mass delivered to constant-concentration receivers (the boundary takes it)
    integer(I4B) :: ritprint = 0 !< E3-ITERATE printed

  contains

    procedure :: allocate_arrays => gwtfmi_allocate_arrays
    procedure :: allocate_gwfpackages => gwtfmi_allocate_gwfpackages
    procedure :: allocate_scalars => gwtfmi_allocate_scalars
    procedure :: deallocate_gwfpackages => gwtfmi_deallocate_gwfpackages
    procedure :: fmi_rp
    procedure :: fmi_ad
    procedure :: fmi_fc
    procedure :: fmi_cq
    procedure :: fmi_bd
    procedure :: fmi_ot_flow
    procedure :: fmi_da => gwtfmi_da
    procedure :: gwfsatold
    procedure :: initialize_gwfterms_from_bfr
    procedure :: initialize_gwfterms_from_gwfbndlist
    procedure :: source_options => gwtfmi_source_options
    procedure :: set_aptbudobj_pointer
    procedure :: source_packagedata => gwtfmi_source_packagedata
    procedure :: set_active_status
    procedure :: route_build
    procedure :: route_fc_reset
    procedure :: route_fc
    procedure :: route_cq_reset
    procedure :: route_cq_faces
    procedure :: route_apt_final
    procedure :: route_cq

  end type TspFmiType

contains

  !> @brief Create a new FMI object
  !<
  subroutine fmi_cr(fmiobj, name_model, input_mempath, inunit, iout, eqnsclfac, &
                    depvartype)
    ! -- dummy
    type(TspFmiType), pointer :: fmiobj
    character(len=*), intent(in) :: name_model
    character(len=*), intent(in) :: input_mempath
    integer(I4B), intent(in) :: inunit
    integer(I4B), intent(in) :: iout
    real(DP), intent(in), pointer :: eqnsclfac !< governing equation scale factor
    character(len=LENVARNAME), intent(in) :: depvartype
    !
    ! -- Create the object
    allocate (fmiobj)
    !
    ! -- create name and memory path
    call fmiobj%set_names(1, name_model, 'FMI', 'FMI', input_mempath)
    fmiobj%text = text
    !
    ! -- Allocate scalars
    call fmiobj%allocate_scalars()
    !
    ! -- Set variables
    fmiobj%inunit = inunit
    fmiobj%iout = iout
    !
    ! -- Assign label based on dependent variable
    fmiobj%depvartype = depvartype
    !
    ! -- Store pointer to governing equation scale factor
    fmiobj%eqnsclfac => eqnsclfac
  end subroutine fmi_cr

  !> @brief Read and prepare
  !<
  subroutine fmi_rp(this, inmvr)
    ! -- modules
    use TdisModule, only: kper, kstp
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: inmvr
    ! -- local
    ! -- formats
    !
    ! --Check to make sure MVT Package is active if mvr flows are available.
    !   This cannot be checked until RP because exchange doesn't set a pointer
    !   to mvrbudobj until exg_ar().
    if (kper * kstp == 1) then
      if (associated(this%mvrbudobj) .and. inmvr == 0) then
        write (errmsg, '(a)') 'GWF water mover is active but the GWT MVT &
          &package has not been specified.  activate GWT MVT package.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
      if (.not. associated(this%mvrbudobj) .and. inmvr > 0) then
        write (errmsg, '(a)') 'GWF water mover terms are not available &
          &but the GWT MVT package has been activated.  Activate GWF-GWT &
          &exchange or specify GWFMOVER in FMI PACKAGEDATA.'
        call store_error(errmsg, terminate=.TRUE.)
      end if
    end if
  end subroutine fmi_rp

  !> @brief Advance routine for FMI object
  !<
  subroutine fmi_ad(this, cnew)
    ! -- modules
    use ConstantsModule, only: DHDRY
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(inout), dimension(:) :: cnew
    ! -- local
    integer(I4B) :: n
    character(len=*), parameter :: fmtdry = &
     &"(/1X,'WARNING: DRY CELL ENCOUNTERED AT ',a,';  RESET AS INACTIVE &
     &WITH DRY CONCENTRATION = ', G13.5)"
    character(len=*), parameter :: fmtrewet = &
     &"(/1X,'DRY CELL REACTIVATED AT ', a,&
     &' WITH STARTING CONCENTRATION =',G13.5)"
    !
    ! -- Set flag to indicated that flows are being updated.  For the case where
    !    flows may be reused (only when flows are read from a file) then set
    !    the flag to zero to indicate that flows were not updated
    this%iflowsupdated = 1
    !
    ! -- If reading flows from a budget file, read the next set of records
    if (this%iubud /= 0) then
      call this%advance_bfr()
    end if
    !
    ! -- If reading heads from a head file, read the next set of records
    if (this%iuhds /= 0) then
      call this%advance_hfr()
    end if
    !
    ! -- If mover flows are being read from file, read the next set of records
    if (this%iumvr /= 0) then
      call this%mvrbudobj%bfr_advance(this%dis, this%iout)
    end if
    !
    ! -- If advanced package flows are being read from file, read the next set of records
    if (this%flows_from_file .and. this%inunit /= 0) then
      do n = 1, size(this%aptbudobj)
        call this%aptbudobj(n)%ptr%bfr_advance(this%dis, this%iout)
      end do
    end if
    !
    ! -- set inactive transport cell status
    if (this%idryinactive /= 0) then
      call this%set_active_status(cnew)
      ! -- E v3: this step's routing table, from this step's flows and active set
      if (this%iroute /= 0) call this%route_build()
    end if
  end subroutine fmi_ad

  !> @brief Calculate coefficients and fill matrix and rhs terms associated
  !! with FMI object
  !<
  subroutine fmi_fc(this, nodes, cold, nja, matrix_sln, idxglo, rhs)
    ! -- dummy
    class(TspFmiType) :: this
    integer, intent(in) :: nodes
    real(DP), intent(in), dimension(nodes) :: cold
    integer(I4B), intent(in) :: nja
    class(MatrixBaseType), pointer :: matrix_sln
    integer(I4B), intent(in), dimension(nja) :: idxglo
    real(DP), intent(inout), dimension(nodes) :: rhs
    ! -- local
    integer(I4B) :: n, idiag, idiag_sln
    real(DP) :: qcorr
    !
    ! -- Calculate the flow imbalance error and make a correction for it
    if (this%iflowerr /= 0) then
      !
      ! -- Correct the transport solution for the flow imbalance by adding
      !    the flow residual to the diagonal
      do n = 1, nodes
        idiag = this%dis%con%ia(n)
        idiag_sln = idxglo(idiag)
        !call matrix_sln%add_value_pos(idiag_sln, -this%gwfflowja(idiag))
        qcorr = -this%gwfflowja(idiag) * this%eqnsclfac
        call matrix_sln%add_value_pos(idiag_sln, qcorr)
      end do
    end if
  end subroutine fmi_fc

  !> @brief Calculate flow correction
  !!
  !! Where there is a flow imbalance for a given cell, a correction may be
  !! applied if selected
  !<
  subroutine fmi_cq(this, cnew, flowja)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(in), dimension(:) :: cnew
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: idiag
    real(DP) :: rate
    !
    ! -- If not adding flow error correction, return
    if (this%iflowerr /= 0) then
      !
      ! -- Accumulate the flow correction term
      do n = 1, this%dis%nodes
        rate = DZERO
        idiag = this%dis%con%ia(n)
        if (this%ibound(n) > 0) then
          rate = -this%gwfflowja(idiag) * cnew(n) * this%eqnsclfac
        end if
        this%flowcorrect(n) = rate
        flowja(idiag) = flowja(idiag) + rate
      end do
    end if
  end subroutine fmi_cq

  !> @brief Calculate budget terms associated with FMI object
  !<
  subroutine fmi_bd(this, isuppress_output, model_budget)
    ! -- modules
    use TdisModule, only: delt
    use BudgetModule, only: BudgetType, rate_accumulator
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: isuppress_output
    type(BudgetType), intent(inout) :: model_budget
    ! -- local
    real(DP) :: rin
    real(DP) :: rout
    !
    ! -- flow correction
    if (this%iflowerr /= 0) then
      call rate_accumulator(this%flowcorrect, rin, rout)
      call model_budget%addentry(rin, rout, delt, budtxt(2), isuppress_output)
    end if
  end subroutine fmi_bd

  !> @brief Save budget terms associated with FMI object
  !<
  subroutine fmi_ot_flow(this, icbcfl, icbcun)
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: icbcfl
    integer(I4B), intent(in) :: icbcun
    ! -- local
    integer(I4B) :: ibinun
    integer(I4B) :: iprint, nvaluesp, nwidthp
    character(len=1) :: cdatafmp = ' ', editdesc = ' '
    real(DP) :: dinact
    !
    ! -- Set unit number for binary output
    if (this%ipakcb < 0) then
      ibinun = icbcun
    elseif (this%ipakcb == 0) then
      ibinun = 0
    else
      ibinun = this%ipakcb
    end if
    if (icbcfl == 0) ibinun = 0
    !
    ! -- Do not save flow corrections if not active
    if (this%iflowerr == 0) ibinun = 0
    !
    ! -- Record the storage rates if requested
    if (ibinun /= 0) then
      iprint = 0
      dinact = DZERO
      !
      ! -- flow correction
      call this%dis%record_array(this%flowcorrect, this%iout, iprint, -ibinun, &
                                 budtxt(2), cdatafmp, nvaluesp, &
                                 nwidthp, editdesc, dinact)
    end if
  end subroutine fmi_ot_flow

  !> @brief Deallocate variables
  !!
  !! Deallocate memory associated with FMI object
  !<
  subroutine gwtfmi_da(this)
    ! -- modules
    use MemoryManagerModule, only: mem_deallocate
    ! -- dummy
    class(TspFmiType) :: this
    ! -- todo: finalize hfr and bfr either here or in a finalize routine
    !
    ! -- deallocate any memory stored with gwfpackages
    call this%deallocate_gwfpackages()
    !
    ! -- deallocate fmi arrays
    if (associated(this%datp)) then
      deallocate (this%datp)
      deallocate (this%gwfpackages)
      deallocate (this%flowpacknamearray)
      call mem_deallocate(this%iatp)
      call mem_deallocate(this%igwfmvrterm)
    end if

    deallocate (this%aptbudobj)
    call mem_deallocate(this%flowcorrect)
    call mem_deallocate(this%ibdgwfsat0)
    if (this%flows_from_file) then
      call mem_deallocate(this%gwfstrgss)
      call mem_deallocate(this%gwfstrgsy)
    end if
    !
    ! -- special treatment, these could be from mem_checkin
    call mem_deallocate(this%gwfhead, 'GWFHEAD', this%memoryPath)
    call mem_deallocate(this%gwfsat, 'GWFSAT', this%memoryPath)
    call mem_deallocate(this%gwfspdis, 'GWFSPDIS', this%memoryPath)
    call mem_deallocate(this%gwfflowja, 'GWFFLOWJA', this%memoryPath)
    !
    ! -- deallocate scalars
    call mem_deallocate(this%flows_from_file)
    call mem_deallocate(this%iflowsupdated)
    call mem_deallocate(this%iflowerr)
    call mem_deallocate(this%igwfstrgss)
    call mem_deallocate(this%igwfstrgsy)
    call mem_deallocate(this%iubud)
    call mem_deallocate(this%iuhds)
    call mem_deallocate(this%iumvr)
    call mem_deallocate(this%iugrb)
    call mem_deallocate(this%nflowpack)
    call mem_deallocate(this%idryinactive)
    !
    ! -- deallocate parent
    call this%NumericalPackageType%da()
  end subroutine gwtfmi_da

  !> @ brief Allocate scalars
  !!
  !! Allocate scalar variables for an FMI object
  !<
  subroutine gwtfmi_allocate_scalars(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate, mem_setptr
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    !
    ! -- allocate scalars in parent
    call this%FlowModelInterfaceType%allocate_scalars()
    !
    ! -- Allocate
    call mem_allocate(this%iflowerr, 'IFLOWERR', this%memoryPath)
    !
    ! -- Although not a scalar, allocate the advanced package transport
    !    budget object to zero so that it can be dynamically resized later
    allocate (this%aptbudobj(0))
    !
    ! -- Initialize
    this%iflowerr = 0
  end subroutine gwtfmi_allocate_scalars

  !> @ brief Allocate arrays for FMI object
  !!
  !!  Method to allocate arrays for the FMI package.
  !<
  subroutine gwtfmi_allocate_arrays(this, nodes)
    use MemoryManagerModule, only: mem_allocate
    ! -- modules
    use ConstantsModule, only: DZERO
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: nodes
    ! -- local
    integer(I4B) :: n
    !
    ! -- allocate parent arrays
    call this%FlowModelInterfaceType%allocate_arrays(nodes)
    !
    ! -- Allocate variables needed for all cases
    if (this%iflowerr == 0) then
      call mem_allocate(this%flowcorrect, 1, 'FLOWCORRECT', this%memoryPath)
    else
      call mem_allocate(this%flowcorrect, nodes, 'FLOWCORRECT', this%memoryPath)
    end if
    do n = 1, size(this%flowcorrect)
      this%flowcorrect(n) = DZERO
    end do
  end subroutine gwtfmi_allocate_arrays

  !> @brief Set gwt transport cell status
  !!
  !! Dry GWF cells are treated differently by GWT and GWE.  Transport does not
  !! occur in deactivated GWF cells; however, GWE still simulates conduction
  !! through dry cells.
  !<
  subroutine set_active_status(this, cnew)
    ! -- modules
    use ConstantsModule, only: DHDRY
    ! -- dummy
    class(TspFmiType) :: this
    real(DP), intent(inout), dimension(:) :: cnew
    ! -- local
    integer(I4B) :: n
    integer(I4B) :: m
    integer(I4B) :: ipos
    real(DP) :: crewet, tflow, flownm
    character(len=15) :: nodestr
    ! -- formats
    character(len=*), parameter :: fmtoutmsg1 = &
      "(1x,'WARNING: DRY CELL ENCOUNTERED AT ', a,'; RESET AS INACTIVE WITH &
      &DRY ', a, '=', G13.5)"
    character(len=*), parameter :: fmtoutmsg2 = &
      &"(1x,'DRY CELL REACTIVATED AT', a, 'WITH STARTING', a, '=', G13.5)"
    !
    ! -- patch E: per-step record of the cells this routine deactivates (read by MST for the DRY-CELL budget)
    if (this%izerosatinact /= 0) then
      if (.not. allocated(this%idrynew)) allocate (this%idrynew(this%dis%nodes))
      this%idrynew(:) = 0
    end if
    do n = 1, this%dis%nodes
      ! -- Calculate the ibound-like array that has 0 if saturation
      !    is zero and 1 otherwise
      if (this%gwfsat(n) > DZERO) then
        this%ibdgwfsat0(n) = 1
      else
        this%ibdgwfsat0(n) = 0
      end if
      !
      ! -- Check if active transport cell is inactive for flow
      if (this%ibound(n) > 0) then
        if (this%gwfhead(n) == DHDRY .or. &
            (this%izerosatinact /= 0 .and. this%gwfsat(n) <= DZERO)) then
          ! -- transport cell should be made inactive (patch E: also at zero flow saturation, when opted in)
          this%ibound(n) = 0
          if (this%izerosatinact /= 0) this%idrynew(n) = 1
          cnew(n) = DHDRY
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtoutmsg1) &
            trim(nodestr), trim(adjustl(this%depvartype)), DHDRY
        end if
      end if
    end do
    !
    ! -- if flow cell is dry, then set gwt%ibound = 0 and conc to dry
    do n = 1, this%dis%nodes
      !
      ! -- Convert dry transport cell to active if flow has rewet
      if (cnew(n) == DHDRY) then
        if (this%gwfhead(n) /= DHDRY .and. &
            (this%izerosatinact == 0 .or. this%gwfsat(n) > DZERO)) then
          !
          ! -- obtain weighted concentration/temperature
          crewet = DZERO
          tflow = DZERO
          do ipos = this%dis%con%ia(n) + 1, this%dis%con%ia(n + 1) - 1
            m = this%dis%con%ja(ipos)
            flownm = this%gwfflowja(ipos)
            if (flownm > 0) then
              if (this%ibound(m) /= 0) then
                crewet = crewet + cnew(m) * flownm ! kluge note: apparently no need to multiply flows by eqnsclfac
                tflow = tflow + this%gwfflowja(ipos) !             since it will divide out below anyway
              end if
            end if
          end do
          if (tflow > DZERO) then
            crewet = crewet / tflow
          else
            crewet = DZERO
          end if
          !
          ! -- cell is now wet
          this%ibound(n) = 1
          cnew(n) = crewet
          call this%dis%noder_to_string(n, nodestr)
          write (this%iout, fmtoutmsg2) &
            trim(nodestr), trim(adjustl(this%depvartype)), crewet
        end if
      end if
    end do
  end subroutine set_active_status

  !> @brief Calculate the previous saturation level
  !!
  !! Calculate the groundwater cell head saturation for the end of
  !! the last time step
  !<
  function gwfsatold(this, n, delt) result(satold)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: n
    real(DP), intent(in) :: delt
    ! -- result
    real(DP) :: satold
    ! -- local
    real(DP) :: vcell
    real(DP) :: vnew
    real(DP) :: vold
    !
    ! -- calculate the value
    vcell = this%dis%area(n) * (this%dis%top(n) - this%dis%bot(n))
    vnew = vcell * this%gwfsat(n)
    vold = vnew
    if (this%igwfstrgss /= 0) vold = vold + this%gwfstrgss(n) * delt
    if (this%igwfstrgsy /= 0) vold = vold + this%gwfstrgsy(n) * delt
    satold = vold / vcell
  end function gwfsatold

  !> @ brief Source input options for package
  !<
  subroutine gwtfmi_source_options(this)
    ! -- modules
    use MemoryManagerExtModule, only: mem_set_value
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    logical(LGP) :: found_ipakcb, found_flowerr
    character(len=*), parameter :: fmtisvflow = &
      "(4x,'CELL-BY-CELL FLOW INFORMATION WILL BE SAVED TO BINARY FILE &
      &WHENEVER ICBCFL IS NOT ZERO AND FLOW IMBALANCE CORRECTION ACTIVE.')"
    character(len=*), parameter :: fmtifc = &
      &"(4x,'MASS WILL BE ADDED OR REMOVED TO COMPENSATE FOR FLOW IMBALANCE.')"

    write (this%iout, '(1x,a)') 'PROCESSING FMI OPTIONS'

    call mem_set_value(this%ipakcb, 'SAVE_FLOWS', this%input_mempath, &
                       found_ipakcb)
    call mem_set_value(this%iflowerr, 'IMBALANCECORRECT', this%input_mempath, &
                       found_flowerr)

    if (found_ipakcb) then
      this%ipakcb = -1
      write (this%iout, fmtisvflow)
    end if
    if (found_flowerr) write (this%iout, fmtifc)

    write (this%iout, '(1x,a)') 'END OF FMI OPTIONS'
  end subroutine gwtfmi_source_options

  !> @ brief Source input options for package
  !<
  subroutine gwtfmi_source_packagedata(this)
    ! -- modules
    use MemoryManagerModule, only: mem_setptr
    use MemoryManagerExtModule, only: mem_set_value, memorystore_release
    use CharacterStringModule, only: CharacterStringType
    use OpenSpecModule, only: ACCESS, FORM
    use ConstantsModule, only: LINELENGTH, DEM6, LENPACKAGENAME
    use InputOutputModule, only: getunit, openfile, urdaux
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: flowtypes
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: fileops
    type(CharacterStringType), dimension(:), contiguous, &
      pointer :: fnames
    type(BudgetObjectType), pointer :: budobjptr
    type(BudObjPtrArray), dimension(:), allocatable :: tmpbudobj
    character(len=LINELENGTH) :: flowtype, fileop, fname
    integer(I4B) :: iapt, inunit, n, i
    logical(LGP) :: exist

    iapt = 0

    call mem_setptr(flowtypes, 'FLOWTYPE', this%input_mempath)
    call mem_setptr(fileops, 'FILEIN', this%input_mempath)
    call mem_setptr(fnames, 'FNAME', this%input_mempath)

    do n = 1, size(flowtypes)
      flowtype = flowtypes(n)
      fileop = fileops(n)
      fname = fnames(n)

      inquire (file=trim(fname), exist=exist)
      if (.not. exist) then
        call store_error('Could not find file '//trim(fname))
        cycle
      end if

      if (fileop /= 'FILEIN') then
        call store_error('Unexpected packagedata input keyword read: "' &
                         //trim(fileop)//'".')
        cycle
      end if

      select case (flowtype)
      case ('GWFBUDGET')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'OLD')
        this%iubud = inunit
        call this%initialize_bfr()
      case ('GWFHEAD')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'OLD')
        this%iuhds = inunit
        call this%initialize_hfr()
      case ('GWFMOVER')
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'OLD')
        this%iumvr = inunit
        call budgetobject_cr_bfr(this%mvrbudobj, 'MVT', this%iumvr, &
                                 this%iout)
        call this%mvrbudobj%fill_from_bfr(this%dis, this%iout)
      case default
        !
        ! --expand the size of aptbudobj, which stores a pointer to the budobj
        allocate (tmpbudobj(iapt))
        do i = 1, size(this%aptbudobj)
          tmpbudobj(i)%ptr => this%aptbudobj(i)%ptr
        end do
        deallocate (this%aptbudobj)
        allocate (this%aptbudobj(iapt + 1))
        do i = 1, size(tmpbudobj)
          this%aptbudobj(i)%ptr => tmpbudobj(i)%ptr
        end do
        deallocate (tmpbudobj)
        !
        ! -- Open the budget file and start filling it
        iapt = iapt + 1
        inunit = getunit()
        call openfile(inunit, this%iout, fname, 'DATA(BINARY)', FORM, &
                      ACCESS, 'OLD')
        call budgetobject_cr_bfr(budobjptr, flowtype, inunit, &
                                 this%iout, colconv2=['GWF             '])
        call budobjptr%fill_from_bfr(this%dis, this%iout)
        this%aptbudobj(iapt)%ptr => budobjptr
      end select
    end do

    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if

    call memorystore_release('FLOWTYPE', this%input_mempath)
    call memorystore_release('FILEIN', this%input_mempath)
    call memorystore_release('FNAME', this%input_mempath)
  end subroutine gwtfmi_source_packagedata

  !> @brief Set the pointer to a budget object
  !!
  !! An advanced transport can pass in a name and a
  !! pointer budget object, and this routine will look through the budget
  !! objects managed by FMI and point to the one with the same name, such as
  !! LAK-1, SFR-1, etc.
  !<
  subroutine set_aptbudobj_pointer(this, name, budobjptr)
    ! -- modules
    class(TspFmiType) :: this
    ! -- dumm
    character(len=*), intent(in) :: name
    type(BudgetObjectType), pointer :: budobjptr
    ! -- local
    integer(I4B) :: i
    !
    ! -- find and set the pointer
    do i = 1, size(this%aptbudobj)
      if (this%aptbudobj(i)%ptr%name == name) then
        budobjptr => this%aptbudobj(i)%ptr
        exit
      end if
    end do
  end subroutine set_aptbudobj_pointer

  !> @brief Initialize the groundwater flow terms based on the budget file
  !! reader
  !!
  !! Initialize terms and figure out how many different terms and packages
  !! are contained within the file
  !<
  subroutine initialize_gwfterms_from_bfr(this)
    ! -- modules
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: nflowpack
    integer(I4B) :: i, ip
    integer(I4B) :: naux
    logical :: found_flowja
    logical :: found_dataspdis
    logical :: found_datasat
    logical :: found_stoss
    logical :: found_stosy
    integer(I4B), dimension(:), allocatable :: imap
    !
    ! -- Calculate the number of gwf flow packages
    allocate (imap(this%bfr%nbudterms))
    imap(:) = 0
    nflowpack = 0
    found_flowja = .false.
    found_dataspdis = .false.
    found_datasat = .false.
    found_stoss = .false.
    found_stosy = .false.
    do i = 1, this%bfr%nbudterms
      select case (trim(adjustl(this%bfr%budtxtarray(i))))
      case ('FLOW-JA-FACE')
        found_flowja = .true.
      case ('DATA-SPDIS')
        found_dataspdis = .true.
      case ('DATA-SAT')
        found_datasat = .true.
      case ('STO-SS')
        found_stoss = .true.
        this%igwfstrgss = 1
      case ('STO-SY')
        found_stosy = .true.
        this%igwfstrgsy = 1
      case default
        nflowpack = nflowpack + 1
        imap(i) = 1
      end select
    end do
    !
    ! -- allocate gwfpackage arrays (gwfpackages, iatp, datp, ...)
    call this%allocate_gwfpackages(nflowpack)
    !
    ! -- Copy the package name and aux names from budget file reader
    !    to the gwfpackages derived-type variable
    ip = 1
    do i = 1, this%bfr%nbudterms
      if (imap(i) == 0) cycle
      call this%gwfpackages(ip)%set_name(this%bfr%dstpackagenamearray(i), &
                                         this%bfr%budtxtarray(i))
      naux = this%bfr%nauxarray(i)
      call this%gwfpackages(ip)%set_auxname(naux, &
                                            this%bfr%auxtxtarray(1:naux, i))
      ip = ip + 1
    end do
    !
    ! -- Copy just the package names for the boundary packages into
    !    the flowpacknamearray
    ip = 1
    do i = 1, size(imap)
      if (imap(i) == 1) then
        this%flowpacknamearray(ip) = this%bfr%dstpackagenamearray(i)
        ip = ip + 1
      end if
    end do
    !
    ! -- Error if specific discharge, saturation or flowja not found
    if (.not. found_dataspdis) then
      write (errmsg, '(a)') 'Specific discharge not found in &
                            &budget file. SAVE_SPECIFIC_DISCHARGE and &
                            &SAVE_FLOWS must be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (.not. found_datasat) then
      write (errmsg, '(a)') 'Saturation not found in &
                            &budget file. SAVE_SATURATION and &
                            &SAVE_FLOWS must be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (.not. found_flowja) then
      write (errmsg, '(a)') 'FLOWJA not found in &
                            &budget file. SAVE_FLOWS must &
                            &be activated in the NPF package.'
      call store_error(errmsg)
    end if
    if (count_errors() > 0) then
      call store_error_filename(this%input_fname)
    end if
  end subroutine initialize_gwfterms_from_bfr

  !> @brief Initialize groundwater flow terms from the groundwater budget
  !!
  !! Flows are coming from a gwf-gwt exchange object
  !<
  subroutine initialize_gwfterms_from_gwfbndlist(this)
    ! -- modules
    use BndModule, only: BndType, GetBndFromList
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: ngwfpack
    integer(I4B) :: ngwfterms
    integer(I4B) :: ip
    integer(I4B) :: imover
    integer(I4B) :: ntomvr
    integer(I4B) :: iterm
    character(len=LENPACKAGENAME) :: budtxt
    class(BndType), pointer :: packobj => null()
    !
    ! -- determine size of gwf terms
    ngwfpack = this%gwfbndlist%Count()
    !
    ! -- Count number of to-mvr terms, but do not include advanced packages
    !    as those mover terms are not losses from the cell, but rather flows
    !    within the advanced package
    ntomvr = 0
    do ip = 1, ngwfpack
      packobj => GetBndFromList(this%gwfbndlist, ip)
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        ntomvr = ntomvr + 1
      end if
    end do
    !
    ! -- Allocate arrays in fmi of size ngwfterms, which is the number of
    !    packages plus the number of packages with mover terms.
    ngwfterms = ngwfpack + ntomvr
    call this%allocate_gwfpackages(ngwfterms)
    !
    ! -- Assign values in the fmi package
    iterm = 1
    do ip = 1, ngwfpack
      !
      ! -- set and store names
      packobj => GetBndFromList(this%gwfbndlist, ip)
      budtxt = adjustl(packobj%text)
      call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
      this%flowpacknamearray(iterm) = packobj%packName
      call this%gwfpackages(iterm)%set_auxname(packobj%naux, &
                                               packobj%auxname)
      iterm = iterm + 1
      !
      ! -- if this package has a mover associated with it, then add another
      !    term that corresponds to the mover flows
      imover = packobj%imover
      if (packobj%isadvpak /= 0) imover = 0
      if (imover /= 0) then
        budtxt = trim(adjustl(packobj%text))//'-TO-MVR'
        call this%gwfpackages(iterm)%set_name(packobj%packName, budtxt)
        this%flowpacknamearray(iterm) = packobj%packName
        call this%gwfpackages(iterm)%set_auxname(packobj%naux, &
                                                 packobj%auxname)
        this%igwfmvrterm(iterm) = 1
        iterm = iterm + 1
      end if
    end do
  end subroutine initialize_gwfterms_from_gwfbndlist

  !> @brief Initialize an array for storing PackageBudget objects.
  !!
  !! This routine allocates gwfpackages (an array of PackageBudget
  !! objects) to the proper size and initializes member variables.
  !<
  subroutine gwtfmi_allocate_gwfpackages(this, ngwfterms)
    ! -- modules
    use ConstantsModule, only: LENMEMPATH
    use MemoryManagerModule, only: mem_allocate
    ! -- dummy
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: ngwfterms
    ! -- local
    integer(I4B) :: n
    character(len=LENMEMPATH) :: memPath
    !
    ! -- direct allocate
    allocate (this%gwfpackages(ngwfterms))
    allocate (this%flowpacknamearray(ngwfterms))
    allocate (this%datp(ngwfterms))
    !
    ! -- mem_allocate
    call mem_allocate(this%iatp, ngwfterms, 'IATP', this%memoryPath)
    call mem_allocate(this%igwfmvrterm, ngwfterms, 'IGWFMVRTERM', this%memoryPath)
    !
    ! -- initialize
    this%nflowpack = ngwfterms
    do n = 1, this%nflowpack
      this%iatp(n) = 0
      this%igwfmvrterm(n) = 0
      this%flowpacknamearray(n) = ''
      !
      ! -- Create a mempath for each individual flow package data set
      !    of the form, MODELNAME/FMI-FTn
      write (memPath, '(a, i0)') trim(this%memoryPath)//'-FT', n
      call this%gwfpackages(n)%initialize(memPath)
    end do
  end subroutine gwtfmi_allocate_gwfpackages

  !> @brief Deallocate memory
  !!
  !! Deallocate memory that stores the gwfpackages array
  !<
  subroutine gwtfmi_deallocate_gwfpackages(this)
    ! -- modules
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: n
    !
    ! -- initialize
    do n = 1, this%nflowpack
      call this%gwfpackages(n)%da()
    end do
  end subroutine gwtfmi_deallocate_gwfpackages

  !> @brief E v3 (rafiei-vahid fork): build this step's routing table
  !!
  !! An E-inactive cell (ibound 0 under DRY_CELL_ZERO_SAT_INACTIVE) holds no water. The water it receives leaves
  !! by its outflows: faces to active cells, faces to other E-inactive cells (a chain, composed top-down), and
  !! sinks (non-advanced GWF packages with outflow at the cell). The mass that arrives with the water is split
  !! over the FINAL destinations by flow share. A cell whose EFFECTIVE outflow (outflows that deliver) is zero is
  !! isolated (rstat 2): its active senders keep their mass (D1), and no chain routes mass into it. A routing
  !! cycle among E-inactive cells, and an advanced package gaining water from an E-inactive cell, are FATAL.
  !<
  subroutine route_build(this)
    ! -- dummy
    class(TspFmiType) :: this
    ! -- local
    integer(I4B) :: nodes, n, m, j, ipos, ip, i, k, c, d, e, nd, nq, iq, nout, nleft
    integer(I4B), dimension(:), allocatable :: oia, otype, otgt, indeg, cursor, kq, topo, opos, hasapt
    real(DP), dimension(:), allocatable :: oq, effout
    real(DP) :: q, f
    type(RouteDestType), dimension(:), allocatable :: fin
    character(len=LINELENGTH) :: cellstr
    !
    nodes = this%dis%nodes
    if (.not. allocated(this%rstat)) then
      allocate (this%rstat(nodes), this%rpos(nodes), this%raptlag(nodes))
      this%raptlag(:) = DZERO
    end if
    this%rstat(:) = 0
    this%rpos(:) = 0
    !
    ! -- 1. the E-inactive cells
    this%nri = 0
    do n = 1, nodes
      if (this%ibound(n) == 0) this%nri = this%nri + 1
    end do
    if (allocated(this%rinode)) deallocate (this%rinode)
    allocate (this%rinode(this%nri))
    k = 0
    do n = 1, nodes
      if (this%ibound(n) == 0) then
        k = k + 1
        this%rinode(k) = n
        this%rpos(n) = k
        this%rstat(n) = 1
      end if
    end do
    !
    ! -- 2. direct outflows of each E-inactive cell: 1 face to an active cell, 3 face to an E-inactive cell,
    !       2 sink (a non-advanced GWF package taking water out at the cell)
    allocate (oia(this%nri + 1))
    oia(:) = 0
    do k = 1, this%nri
      m = this%rinode(k)
      do ipos = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        if (this%gwfflowja(ipos) < DZERO) oia(k + 1) = oia(k + 1) + 1
      end do
    end do
    do ip = 1, this%nflowpack
      do i = 1, this%gwfpackages(ip)%nbound
        n = this%gwfpackages(ip)%nodelist(i)
        if (n <= 0) cycle
        if (this%rpos(n) == 0) cycle
        q = this%gwfpackages(ip)%get_flow(i)
        if (q < DZERO) then
          if (this%iatp(ip) /= 0) then
            call this%dis%noder_to_string(n, cellstr)
            write (errmsg, '(a,a,a,a,a,es15.6,a)') &
              'DRY_CELL_ROUTE_MASS: the advanced package ', trim(this%flowpacknamearray(ip)), &
              ' GAINS water from the E-inactive cell ', trim(cellstr), ' (', q, &
              '); routing mass into a feature is not implemented.'
            call store_error(errmsg, terminate=.TRUE.)
          end if
          oia(this%rpos(n) + 1) = oia(this%rpos(n) + 1) + 1
        end if
      end do
    end do
    oia(1) = 1
    do k = 1, this%nri
      oia(k + 1) = oia(k + 1) + oia(k)
    end do
    nout = oia(this%nri + 1) - 1
    allocate (otype(max(nout, 1)), otgt(max(nout, 1)), oq(max(nout, 1)), opos(max(nout, 1)))
    allocate (indeg(max(this%nri, 1)), cursor(max(this%nri, 1)), hasapt(max(this%nri, 1)))
    indeg(:) = 0
    opos(:) = 0
    hasapt(:) = 0
    if (this%nri > 0) cursor(1:this%nri) = oia(1:this%nri)
    do k = 1, this%nri
      m = this%rinode(k)
      do ipos = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        q = this%gwfflowja(ipos)
        if (q < DZERO) then
          j = this%dis%con%ja(ipos)
          e = cursor(k)
          cursor(k) = e + 1
          oq(e) = -q
          if (this%rpos(j) /= 0) then
            otype(e) = 3
            otgt(e) = this%rpos(j)
            indeg(this%rpos(j)) = indeg(this%rpos(j)) + 1
          else
            otype(e) = 1
            otgt(e) = j
            opos(e) = this%dis%con%isym(ipos)
          end if
        end if
      end do
    end do
    do ip = 1, this%nflowpack
      do i = 1, this%gwfpackages(ip)%nbound
        n = this%gwfpackages(ip)%nodelist(i)
        if (n <= 0) cycle
        if (this%rpos(n) == 0) cycle
        q = this%gwfpackages(ip)%get_flow(i)
        if (q < DZERO) then
          k = this%rpos(n)
          e = cursor(k)
          cursor(k) = e + 1
          otype(e) = 2
          otgt(e) = n
          oq(e) = -q
        else if (q > DZERO .and. this%iatp(ip) /= 0) then
          hasapt(this%rpos(n)) = 1
        end if
      end do
    end do
    !
    ! -- 3. upstream-first order over E-inactive -> E-inactive faces (Kahn); a cycle is FATAL
    allocate (kq(max(this%nri, 1)), topo(max(this%nri, 1)))
    nq = 0
    do k = 1, this%nri
      if (indeg(k) == 0) then
        nq = nq + 1
        kq(nq) = k
      end if
    end do
    iq = 1
    c = 0
    do while (iq <= nq)
      k = kq(iq)
      iq = iq + 1
      c = c + 1
      topo(c) = k
      do e = oia(k), oia(k + 1) - 1
        if (otype(e) /= 3) cycle
        indeg(otgt(e)) = indeg(otgt(e)) - 1
        if (indeg(otgt(e)) == 0) then
          nq = nq + 1
          kq(nq) = otgt(e)
        end if
      end do
    end do
    if (c < this%nri) then
      nleft = 0
      do k = 1, this%nri
        if (indeg(k) > 0) then
          nleft = nleft + 1
          if (nleft <= 20) then
            call this%dis%noder_to_string(this%rinode(k), cellstr)
            write (this%iout, '(1x,a,a)') 'E3-CYCLE cell ', trim(cellstr)
          end if
        end if
      end do
      write (errmsg, '(a,i0,a)') 'DRY_CELL_ROUTE_MASS: a routing CYCLE among E-inactive cells (', nleft, &
        ' cells in or below it; the first 20 are printed as E3-CYCLE lines in the listing).'
      call store_error(errmsg, terminate=.TRUE.)
    end if
    !
    ! -- 4. downstream first: effective outflow, isolation, and the final destinations with their fractions
    allocate (effout(max(this%nri, 1)), fin(max(this%nri, 1)))
    do c = this%nri, 1, -1
      k = topo(c)
      m = this%rinode(k)
      effout(k) = DZERO
      do e = oia(k), oia(k + 1) - 1
        if (otype(e) == 3) then
          if (this%rstat(this%rinode(otgt(e))) == 1) effout(k) = effout(k) + oq(e)
        else
          effout(k) = effout(k) + oq(e)
        end if
      end do
      if (effout(k) <= DZERO) then
        this%rstat(m) = 2
        allocate (fin(k)%t(0), fin(k)%n(0), fin(k)%f(0), fin(k)%p(0))
        cycle
      end if
      nd = 0
      do e = oia(k), oia(k + 1) - 1
        if (otype(e) == 3) then
          if (this%rstat(this%rinode(otgt(e))) == 1) nd = nd + size(fin(otgt(e))%t)
        else
          nd = nd + 1
        end if
      end do
      allocate (fin(k)%t(nd), fin(k)%n(nd), fin(k)%f(nd), fin(k)%p(nd))
      d = 0
      do e = oia(k), oia(k + 1) - 1
        f = oq(e) / effout(k)
        if (otype(e) == 3) then
          j = otgt(e)
          if (this%rstat(this%rinode(j)) == 1) then
            do i = 1, size(fin(j)%t)
              d = d + 1
              fin(k)%t(d) = fin(j)%t(i)
              fin(k)%n(d) = fin(j)%n(i)
              fin(k)%f(d) = f * fin(j)%f(i)
              fin(k)%p(d) = fin(j)%p(i)
            end do
          end if
        else
          d = d + 1
          fin(k)%t(d) = otype(e)
          fin(k)%n(d) = otgt(e)
          fin(k)%f(d) = f
          fin(k)%p(d) = opos(e)
        end if
      end do
      ! -- advanced-package mass is delivered after the packages' flow calculation, too late for a constant-
      !    concentration receiver's own balance: such a route stops the run instead of mis-booking
      if (hasapt(k) /= 0) then
        do d = 1, nd
          if (fin(k)%t(d) == 1) then
            if (this%ibound(fin(k)%n(d)) < 0) then
              call this%dis%noder_to_string(fin(k)%n(d), cellstr)
              write (errmsg, '(a,a,a)') 'DRY_CELL_ROUTE_MASS: advanced-package mass routed through an E-inactive cell &
                &would reach the constant-concentration cell ', trim(cellstr), '; not implemented.'
              call store_error(errmsg, terminate=.TRUE.)
            end if
          end if
        end do
      end if
    end do
    !
    ! -- 5. flatten into a CSR over the rinode positions
    if (allocated(this%rdia)) deallocate (this%rdia)
    if (allocated(this%rdtype)) deallocate (this%rdtype, this%rdnode, this%rdfrac, this%rdpos)
    allocate (this%rdia(this%nri + 1))
    this%rdia(1) = 1
    do k = 1, this%nri
      this%rdia(k + 1) = this%rdia(k) + size(fin(k)%t)
    end do
    nd = this%rdia(this%nri + 1) - 1
    allocate (this%rdtype(max(nd, 1)), this%rdnode(max(nd, 1)), this%rdfrac(max(nd, 1)), this%rdpos(max(nd, 1)))
    do k = 1, this%nri
      d = this%rdia(k)
      do i = 1, size(fin(k)%t)
        this%rdtype(d) = fin(k)%t(i)
        this%rdnode(d) = fin(k)%n(i)
        this%rdfrac(d) = fin(k)%f(i)
        this%rdpos(d) = fin(k)%p(i)
        d = d + 1
      end do
    end do
    deallocate (oia, otype, otgt, oq, opos, indeg, cursor, hasapt, kq, topo, effout, fin)
  end subroutine route_build

  !> @brief E v3: clear the advanced-package routing registrations before an outer iteration's fill
  !<
  subroutine route_fc_reset(this)
    class(TspFmiType) :: this
    if (allocated(this%raptlag)) this%raptlag(:) = DZERO
  end subroutine route_fc_reset

  !> @brief E v3: deliver the routed mass to the active receivers (explicit, Picard-lagged)
  !!
  !! Called after every package's fill in the same outer iteration: each routable E-inactive cell's inflowing
  !! mass (faces from active cells at the current iterate, first-order upwind; advanced packages as registered
  !! by their fill) is added as a source to its final receivers by flow share. The senders' loss is implicit
  !! in the advection fill, so the two agree at convergence.
  !<
  subroutine route_fc(this, cnew, rhs)
    class(TspFmiType) :: this
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), intent(inout) :: rhs
    integer(I4B) :: k, m, ipos, j, d
    real(DP) :: q, fin, sumd
    !
    sumd = DZERO
    if (.not. allocated(this%rstat)) return
    do k = 1, this%nri
      m = this%rinode(k)
      if (this%rstat(m) /= 1) cycle
      fin = this%raptlag(m)
      do ipos = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        j = this%dis%con%ja(ipos)
        if (this%ibound(j) == 0) cycle
        q = this%gwfflowja(ipos)
        if (q > DZERO) fin = fin + q * this%eqnsclfac * cnew(j)
      end do
      do d = this%rdia(k), this%rdia(k + 1) - 1
        if (this%rdtype(d) == 1) then
          rhs(this%rdnode(d)) = rhs(this%rdnode(d)) - this%rdfrac(d) * fin
          sumd = sumd + this%rdfrac(d) * fin
        end if
      end do
    end do
    this%rsumd_fc = sumd
  end subroutine route_fc

  !> @brief E v3: the face part of the routed mass in flowja, BEFORE the packages' flow calculation
  !!
  !! A constant-concentration package balances the FACE flows of its cell (gwt-cnc.f90 cnc_cq sums the off-diagonal
  !! flowja), so the routed mass must sit at faces, and before that package runs: the sender's removal at the face
  !! (sender, E-inactive cell), and each face-origin delivery at the face (receiver, adjacent E-inactive cell). The
  !! advanced-package deliveries follow in route_cq, at faces too (a constant-concentration receiver of those is FATAL
  !! in route_build).
  !<
  subroutine route_cq_faces(this, cnew, flowja)
    class(TspFmiType) :: this
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B) :: k, m, ipos, j, d, jpos
    real(DP) :: q, mass, fin
    !
    if (.not. allocated(this%rstat)) return
    do k = 1, this%nri
      m = this%rinode(k)
      if (this%rstat(m) /= 1) cycle
      fin = DZERO
      do ipos = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        j = this%dis%con%ja(ipos)
        if (this%ibound(j) == 0) cycle
        q = this%gwfflowja(ipos)
        if (q <= DZERO) cycle
        mass = q * this%eqnsclfac * cnew(j)
        fin = fin + mass
        jpos = this%dis%con%isym(ipos)
        flowja(jpos) = flowja(jpos) - mass
      end do
      do d = this%rdia(k), this%rdia(k + 1) - 1
        if (this%rdtype(d) == 1) flowja(this%rdpos(d)) = flowja(this%rdpos(d)) + this%rdfrac(d) * fin
      end do
    end do
  end subroutine route_cq_faces

  !> @brief E v3: clear the converged advanced-package registrations before a step's flow calculation
  !<
  subroutine route_cq_reset(this)
    class(TspFmiType) :: this
    this%naptc = 0
  end subroutine route_cq_reset

  !> @brief E v3: register one losing advanced-package connection into an E-inactive cell (converged)
  !<
  subroutine route_apt_final(this, ifeat, node, qwater, qmass)
    class(TspFmiType) :: this
    integer(I4B), intent(in) :: ifeat
    integer(I4B), intent(in) :: node
    real(DP), intent(in) :: qwater
    real(DP), intent(in) :: qmass
    if (.not. allocated(this%raptf)) then
      allocate (this%raptf(64), this%raptm(64), this%raptw(64), this%raptq(64))
    end if
    if (this%naptc + 1 > size(this%raptf)) then
      call route_grow_i(this%raptf)
      call route_grow_i(this%raptm)
      call route_grow_d(this%raptw)
      call route_grow_d(this%raptq)
    end if
    this%naptc = this%naptc + 1
    this%raptf(this%naptc) = ifeat
    this%raptm(this%naptc) = node
    this%raptw(this%naptc) = qwater
    this%raptq(this%naptc) = qmass
  end subroutine route_apt_final

  !> @brief E v3: book the routed mass from the CONVERGED iterate, and print the E3 lines
  !!
  !! ROUTED: - at each active sender (the face's q x C_sender), + at each receiver (its share; id2 = the feature
  !! for an advanced-package origin, else 0). ROUTED-SINK: - at the E-inactive cell of each sink. WITHHELD-ISOL:
  !! + at the sender for a face into an isolated cell (kept there, D1), + at the cell for an advanced-package
  !! connection into an isolated cell (the feature loses it, as E v2 does). Each ROUTED entry is also added to
  !! the cell's flowja diagonal, so the per-cell residual stays closed.
  !<
  subroutine route_cq(this, cnew, flowja)
    use TdisModule, only: kper, kstp
    class(TspFmiType) :: this
    real(DP), dimension(:), intent(in) :: cnew
    real(DP), dimension(:), contiguous, intent(inout) :: flowja
    integer(I4B) :: k, m, ipos, j, d, c, niso
    integer(I4B), dimension(:), allocatable :: isoin
    real(DP) :: q, mass, fin, share, sumd, ftot, atot, wsum, qfaces, qapt
    !
    this%nro = 0
    this%nrs = 0
    this%nrw = 0
    this%rcncin = DZERO
    this%rcncout = DZERO
    sumd = DZERO
    ftot = DZERO
    atot = DZERO
    wsum = DZERO
    qfaces = DZERO
    qapt = DZERO
    if (.not. allocated(this%rstat)) return
    allocate (isoin(max(this%nri, 1)))
    isoin(:) = 0
    do k = 1, this%nri
      m = this%rinode(k)
      fin = DZERO
      do ipos = this%dis%con%ia(m) + 1, this%dis%con%ia(m + 1) - 1
        if (this%dis%con%mask(ipos) == 0) cycle
        j = this%dis%con%ja(ipos)
        if (this%ibound(j) == 0) cycle
        q = this%gwfflowja(ipos)
        if (q <= DZERO) cycle
        mass = q * this%eqnsclfac * cnew(j)
        qfaces = qfaces + q
        if (this%rstat(m) == 1) then
          call route_push(this%nro, this%ronode, this%roid2, this%roq, j, 0, -mass)
          if (this%ibound(j) < 0) this%rcncin = this%rcncin + mass
          ftot = ftot + mass
          fin = fin + mass
        else
          call route_push(this%nrw, this%rwnode, this%rwid2, this%rwq, j, 0, mass)
          wsum = wsum + mass
          isoin(k) = 1
        end if
      end do
      if (this%rstat(m) == 1) then
        do d = this%rdia(k), this%rdia(k + 1) - 1
          share = this%rdfrac(d) * fin
          if (this%rdtype(d) == 1) then
            call route_push(this%nro, this%ronode, this%roid2, this%roq, this%rdnode(d), 0, share)
            if (this%ibound(this%rdnode(d)) < 0) this%rcncout = this%rcncout + share
            sumd = sumd + share
          else
            call route_push(this%nrs, this%rsnode, this%rsid2, this%rsq, this%rdnode(d), 0, -share)
          end if
        end do
      end if
    end do
    do c = 1, this%naptc
      m = this%raptm(c)
      k = this%rpos(m)
      mass = this%raptq(c)
      qapt = qapt + this%raptw(c)
      atot = atot + mass
      if (k > 0) then
        if (this%rstat(m) == 1) then
          do d = this%rdia(k), this%rdia(k + 1) - 1
            share = this%rdfrac(d) * mass
            if (this%rdtype(d) == 1) then
              call route_push(this%nro, this%ronode, this%roid2, this%roq, this%rdnode(d), &
                              this%raptf(c), share)
              if (this%ibound(this%rdnode(d)) < 0) this%rcncout = this%rcncout + share
              flowja(this%rdpos(d)) = flowja(this%rdpos(d)) + share
              sumd = sumd + share
            else
              call route_push(this%nrs, this%rsnode, this%rsid2, this%rsq, this%rdnode(d), &
                              this%raptf(c), -share)
            end if
          end do
          cycle
        end if
        isoin(k) = 1
      end if
      call route_push(this%nrw, this%rwnode, this%rwid2, this%rwq, m, this%raptf(c), mass)
      wsum = wsum + mass
    end do
    niso = 0
    do k = 1, this%nri
      if (isoin(k) == 1) niso = niso + 1
    end do
    deallocate (isoin)
    if (this%ritprint == 0) then
      write (this%iout, '(/1x,a)') 'E3-ITERATE converged'
      write (this%iout, '(1x,a)') 'DRY_CELL_ROUTE_MASS: ROUTED, ROUTED-SINK and WITHHELD-ISOL are booked from the CONVERGED &
        &iterate; E3-LAG dROUTED_IN = |deliveries at the converged iterate - deliveries the last outer iteration &
        &applied|; ROUTED_OUT = face mass removed at senders + advanced-package mass into E-inactive cells.'
      this%ritprint = 1
    end if
    write (this%iout, '(1x,a,2(1x,i0),1x,es24.16e3,1x,i0)') 'E3-ISOLATED', kper, kstp, wsum, niso
    write (this%iout, '(1x,a,2(1x,i0),2(1x,es24.16e3))') 'E3-LAG', kper, kstp, abs(sumd - this%rsumd_fc), &
      ftot + atot
    write (this%iout, '(1x,a,2(1x,i0),2(1x,es24.16e3))') 'E3-QROUTED', kper, kstp, qfaces, qapt
  end subroutine route_cq

  !> @brief E v3: append one entry to a routing output list
  !<
  subroutine route_push(nlist, anode, aid2, aq, node, id2, q)
    integer(I4B), intent(inout) :: nlist
    integer(I4B), dimension(:), allocatable, intent(inout) :: anode
    integer(I4B), dimension(:), allocatable, intent(inout) :: aid2
    real(DP), dimension(:), allocatable, intent(inout) :: aq
    integer(I4B), intent(in) :: node
    integer(I4B), intent(in) :: id2
    real(DP), intent(in) :: q
    if (.not. allocated(anode)) then
      allocate (anode(256), aid2(256), aq(256))
    end if
    if (nlist + 1 > size(anode)) then
      call route_grow_i(anode)
      call route_grow_i(aid2)
      call route_grow_d(aq)
    end if
    nlist = nlist + 1
    anode(nlist) = node
    aid2(nlist) = id2
    aq(nlist) = q
  end subroutine route_push

  subroutine route_grow_i(a)
    integer(I4B), dimension(:), allocatable, intent(inout) :: a
    integer(I4B), dimension(:), allocatable :: b
    allocate (b(2 * max(size(a), 32)))
    b(1:size(a)) = a
    call move_alloc(b, a)
  end subroutine route_grow_i

  subroutine route_grow_d(a)
    real(DP), dimension(:), allocatable, intent(inout) :: a
    real(DP), dimension(:), allocatable :: b
    allocate (b(2 * max(size(a), 32)))
    b(1:size(a)) = a
    call move_alloc(b, a)
  end subroutine route_grow_d

end module TspFmiModule
