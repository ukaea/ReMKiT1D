!-----------------------------------------------------------------------------------------------------------------------------------
! This file is part of ReMKiT1D.
!
! ReMKiT1D is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as 
! published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
! ReMKiT1D is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of 
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with ReMKiT1D. If not, see <https://www.gnu.org/licenses/>. 
!
! Copyright 2023 United Kingdom Atomic Energy Authority (stefan.mijin@ukaea.uk)
!-----------------------------------------------------------------------------------------------------------------------------------
submodule (calculation_tree_class) calculation_tree_procedures
!! author: Stefan Mijin 
!! 
!!  Contains module procedures associated with calculation tree and node classes

implicit none

!-----------------------------------------------------------------------------------------------------------------------------------
contains
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initNode(this,additiveMode,constant,leafVarIndex,unaryRealParams,&
                            unaryIntParams,unaryLogicalParams,unaryTransformTag)
    !! Calculation node initialization routine

    class(CalculationNode)                         ,intent(inout)  :: this
    logical                              ,optional ,intent(in) :: additiveMode
    real(rk)                             ,optional ,intent(in) :: constant
    integer(ik)                          ,optional ,intent(in) :: leafVarIndex
    real(rk)    ,dimension(:)            ,optional ,intent(in) :: unaryRealParams
    integer(ik) ,dimension(:)            ,optional ,intent(in) :: unaryIntParams
    logical     ,dimension(:)            ,optional ,intent(in) :: unaryLogicalParams
    character(*)                         ,optional ,intent(in) :: unaryTransformTag

    procedure(realArrayFunctionGenParam) ,pointer  :: unaryTransform

    this%kernel%additiveMode = .false. 

    if (present(additiveMode)) then
        this%kernel%additiveMode = additiveMode
    end if
    this%kernel%constant = real(1,kind=rk) 
    if (this%kernel%additiveMode) this%kernel%constant = 0
    if (present(constant)) this%kernel%constant = constant

    this%kernel%leafVarIndex = 0
    if (present(leafVarIndex)) this%kernel%leafVarIndex = leafVarIndex

    if (present(unaryTransformTag)) then 
        this%kernel%unaryTransformationTag = unaryTransformTag
        if (allocated(this%kernel%unaryTransformationTag)) then 
            call associateFunctionPointer(unaryTransformTag,this%unaryTransform)
            if (present(unaryIntParams)) this%kernel%unaryIntParams = unaryIntParams
            if (present(unaryRealParams)) this%kernel%unaryRealParams = unaryRealParams
            if (present(unaryLogicalParams)) this%kernel%unaryLogicalParams = unaryLogicalParams
        end if
    end if 

    call this%makeDefined()

end subroutine initNode
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initTree(this,additiveMode,constant,leafVarIndex,unaryRealParams,&
        unaryIntParams,unaryLogicalParams,unaryTransformTag)
    !! Calculation tree initialization routine

    class(CalculationTree)                         ,intent(inout)  :: this
    logical                              ,optional ,intent(in) :: additiveMode
    real(rk)                             ,optional ,intent(in) :: constant
    integer(ik)                          ,optional ,intent(in) :: leafVarIndex
    real(rk)    ,dimension(:)            ,optional ,intent(in) :: unaryRealParams
    integer(ik) ,dimension(:)            ,optional ,intent(in) :: unaryIntParams
    logical     ,dimension(:)            ,optional ,intent(in) :: unaryLogicalParams
    character(*)                         ,optional ,intent(in) :: unaryTransformTag

    allocate(this%root)

    call this%root%init(additiveMode,constant,leafVarIndex,unaryRealParams,&
                        unaryIntParams,unaryLogicalParams,unaryTransformTag)

    call this%makeDefined()

end subroutine initTree
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine addChild(this,additiveMode,constant,leafVarIndex,unaryRealParams,&
    unaryIntParams,unaryLogicalParams,unaryTransformTag)
    !! Initialize a child node of this node with given properties

    class(CalculationNode)               ,target   ,intent(inout)  :: this
    logical                              ,optional ,intent(in) :: additiveMode
    real(rk)                             ,optional ,intent(in) :: constant
    integer(ik)                          ,optional ,intent(in) :: leafVarIndex
    real(rk)    ,dimension(:)            ,optional ,intent(in) :: unaryRealParams
    integer(ik) ,dimension(:)            ,optional ,intent(in) :: unaryIntParams
    logical     ,dimension(:)            ,optional ,intent(in) :: unaryLogicalParams
    character(*)                         ,optional ,intent(in) :: unaryTransformTag

    type(CalculationNode) ,pointer :: nodePointer 

    if (.not. associated(this%leftChild)) then
        allocate(this%leftChild)
        call this%leftChild%init(additiveMode,constant,leafVarIndex,unaryRealParams,&
                                unaryIntParams,unaryLogicalParams,unaryTransformTag)
        this%leftChild%parent => this
    else
        nodePointer => this%leftChild
        do 
            if (associated(nodePointer%rightSibling)) then 
                nodePointer => nodePointer%rightSibling
                cycle
            end if
            allocate(nodePointer%rightSibling)
            call nodePointer%rightSibling%init(additiveMode,constant,leafVarIndex,unaryRealParams,&
                                                unaryIntParams,unaryLogicalParams,unaryTransformTag)
            nodePointer%rightSibling%parent => this
            exit
        end do
    end if

end subroutine addChild
!-----------------------------------------------------------------------------------------------------------------------------------
pure module recursive function evaluateNode(this,inputArray) result(res)
    !! Recursively evaluate nodes, using the inputArray variables for leaf values

    class(CalculationNode)        ,intent(in) :: this
    type(RealArray) ,dimension(:) ,intent(in) :: inputArray
    real(rk) ,allocatable ,dimension(:)       :: res

    if (associated(this%leftChild)) then
        res = this%leftChild%evaluate(inputArray)
    else
        !No check here to make sure that the leaf has a valid variable index
        res = inputArray(this%kernel%leafVarIndex)%entry
    end if

    if (this%kernel%additiveMode) then 
        res = res + this%kernel%constant
    else
        res = res * this%kernel%constant
    end if

    if (associated(this%unaryTransform)) &
        res = this%unaryTransform(res,this%kernel%unaryRealParams,this%kernel%unaryIntParams,this%kernel%unaryLogicalParams)

    if (associated(this%rightSibling)) then 
        if (this%parent%kernel%additiveMode) then
            res = res + this%rightSibling%evaluate(inputArray)
        else
            res = res * this%rightSibling%evaluate(inputArray)
        end if
    end if

end function evaluateNode
!-----------------------------------------------------------------------------------------------------------------------------------
pure module function evaluateTree(this,inputArray) result(res)
    !! Call tree's root node evaluate

    class(CalculationTree)        ,intent(in) :: this
    type(RealArray) ,dimension(:) ,intent(in) :: inputArray
    real(rk) ,allocatable ,dimension(:)       :: res

    res = this%root%evaluate(inputArray)

end function evaluateTree
!-----------------------------------------------------------------------------------------------------------------------------------
module function flattenTree(this) result(res)
    !! Flatten tree into FlatTree object

    class(CalculationTree)        ,intent(in) :: this
    type(FlatTree)                            :: res
    type(CalculationNode) ,pointer :: nodePointer 

    integer(ik) :: numNodes ,parentIndex, currentIndex ,i

    if (assertions) call assert(this%isDefined(),"Attempted to flatten undefined CalculationTree")

    numNodes = 1
    nodePointer => this%root

    !Count numNodes
    do 
        if (associated(nodePointer%leftChild)) then
            numNodes = numNodes + 1
            nodePointer => nodePointer%leftChild
            cycle
        end if

        if (associated(nodePointer%rightSibling)) then
            numNodes = numNodes + 1
            nodePointer => nodePointer%rightSibling
            cycle
        end if

        do 
            if (.not. associated(nodePointer%parent)) exit
            if (associated(nodePointer%parent%rightSibling)) then
                nodePointer => nodePointer%parent%rightSibling
                numNodes = numNodes + 1
                exit

            else
                nodePointer => nodePointer%parent
            end if
        end do
        if (.not. associated(nodePointer%parent)) exit
    end do

    allocate(res%kernels(numNodes))
    allocate(res%children(numNodes))
    do i = 1,numNodes
        allocate(res%children(i)%entry(0))
    end do
    allocate(res%parent(numNodes))

    parentIndex = 0
    currentIndex = 1 

    res%kernels(currentIndex) = this%root%kernel
    res%parent(currentIndex) = 0
    nodePointer => this%root
    do 
        if (associated(nodePointer%leftChild)) then
            res%children(currentIndex)%entry = [res%children(currentIndex)%entry,currentIndex+1]
            parentIndex = currentIndex
            currentIndex = currentIndex + 1
            nodePointer => nodePointer%leftChild
            res%kernels(currentIndex) = nodePointer%kernel
            res%parent(currentIndex) = parentIndex
            cycle
        end if

        if (associated(nodePointer%rightSibling)) then
            currentIndex = currentIndex + 1
            res%children(parentIndex)%entry = [res%children(parentIndex)%entry,currentIndex]
            nodePointer => nodePointer%rightSibling
            res%kernels(currentIndex) = nodePointer%kernel
            res%parent(currentIndex) = parentIndex
            cycle
        end if

        do 
            if (.not. associated(nodePointer%parent)) exit
            if (associated(nodePointer%parent%rightSibling)) then
                nodePointer => nodePointer%parent%rightSibling
                parentIndex = res%parent(parentIndex)
                currentIndex = currentIndex + 1
                res%children(parentIndex)%entry = [res%children(parentIndex)%entry,currentIndex]
                res%kernels(currentIndex) = nodePointer%kernel
                res%parent(currentIndex) = parentIndex
                exit

            else
                nodePointer => nodePointer%parent
                parentIndex = res%parent(parentIndex)
            end if
        end do
        if (.not. associated(nodePointer%parent)) exit
    end do
end function flattenTree
!-----------------------------------------------------------------------------------------------------------------------------------
module subroutine initFromFlatTree(this,fTree)
    !! Calculation tree initialization routine using a FlatTree object

    class(CalculationTree)           ,intent(inout)  :: this
    type(FlatTree)                   ,intent(in)     :: fTree

    type(CalculationNode) ,pointer :: nodePointer 

    integer(ik) :: i ,currentIndex ,parentIndex ,siblingIndex

    integer(ik) ,dimension(:) ,allocatable :: indexLookup

    if (assertions) &
    call assert(.not. this%isDefined(),"Cannot initialize CalculationTree from FlatTree if it is already defined")

    call this%init(fTree%kernels(1)%additiveMode,&
                    fTree%kernels(1)%constant,&
                    fTree%kernels(1)%leafVarIndex,&
                    fTree%kernels(1)%unaryRealParams,&
                    fTree%kernels(1)%unaryIntParams,&
                    fTree%kernels(1)%unaryLogicalParams,&
                    fTree%kernels(1)%unaryTransformationTag)

    parentIndex = 0
    currentIndex = 1
    siblingIndex = 1
    nodePointer => this%root
    do 
        do i = 1,size(fTree%children(currentIndex)%entry)
            call nodePointer%addChild(fTree%kernels(fTree%children(currentIndex)%entry(i))%additiveMode,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%constant,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%leafVarIndex,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%unaryRealParams,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%unaryIntParams,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%unaryLogicalParams,&
                                        fTree%kernels(fTree%children(currentIndex)%entry(i))%unaryTransformationTag)
        end do

        if (associated(nodePointer%leftChild)) then 
            nodePointer => nodePointer%leftChild
            parentIndex = currentIndex
            currentIndex = fTree%children(currentIndex)%entry(1)
            cycle
        end if

        if (associated(nodePointer%rightSibling)) then 
            nodePointer => nodePointer%rightSibling
            indexLookup = findIndices(fTree%children(parentIndex)%entry == currentIndex)
            siblingIndex = indexLookup(1) + 1
            currentIndex = fTree%children(parentIndex)%entry(siblingIndex)
            cycle
        end if

        do 
            if (.not. associated(nodePointer%parent)) exit
            if (associated(nodePointer%parent%rightSibling)) then
                nodePointer => nodePointer%parent%rightSibling
                currentIndex = parentIndex
                parentIndex = fTree%parent(parentIndex)
                indexLookup = findIndices(fTree%children(parentIndex)%entry == currentIndex)
                siblingIndex = indexLookup(1) + 1
                currentIndex =  fTree%children(parentIndex)%entry(siblingIndex)
                exit

            else
                nodePointer => nodePointer%parent
                currentIndex = parentIndex
                parentIndex = fTree%parent(parentIndex)
            end if
        end do
        if (.not. associated(nodePointer%parent)) exit

    end do

end subroutine initFromFlatTree
!-----------------------------------------------------------------------------------------------------------------------------------
pure module recursive subroutine destroyNode(this)

    class(CalculationNode),intent(inout) :: this

    if (associated(this%leftChild)) call this%leftChild%destroy()
    if (associated(this%rightSibling)) call this%rightSibling%destroy()
    nullify(this%leftChild)
    nullify(this%rightSibling)
    nullify(this%parent)
    nullify(this%unaryTransform)

end subroutine destroyNode
!-----------------------------------------------------------------------------------------------------------------------------------
elemental module subroutine finalizeCalculationTree(this) 

    type(CalculationTree) ,intent(inout) :: this

    call this%root%destroy()
    nullify(this%root)

end subroutine finalizeCalculationTree 
!-----------------------------------------------------------------------------------------------------------------------------------
end submodule calculation_tree_procedures
!-----------------------------------------------------------------------------------------------------------------------------------
