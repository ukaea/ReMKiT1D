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
module calculation_tree_class
    !! author: Stefan Mijin 
    !! 
    !! Houses the calculation tree class and the relevant node class. 

    use data_kinds                  ,only: rk, ik
    use runtime_constants           ,only: debugging, assertions
    use god_objects                 ,only: Object
    use assertion_utility           ,only: assert, assertIdentical, assertPure
    use basic_interfaces            ,only: realArrayFunctionGenParam
    use support_functions
    use unary_transforms
    use support_types

    implicit none

    type ,public :: CalculationKernel
        !! Kernel containing calculation node properties
        logical :: additiveMode !! If true will use additive mode for reducing the results of its children. Defaults to false. 

        real(rk) :: constant !! Constant component, defaults to 0 in additive mode and 1 in multiplicative.

        integer(ik) :: leafVarIndex !! Index of variable to be used a the child's result if this node is a leaf

        real(rk)    ,allocatable ,dimension(:) :: unaryRealParams !! Optional real parameters of unary transformation
        integer(ik) ,allocatable ,dimension(:) :: unaryIntParams !! Optional integer parameters of unary transformation
        logical     ,allocatable ,dimension(:) :: unaryLogicalParams !! Optional logical parameters of unary transformation

        character(:) ,allocatable :: unaryTransformationTag 
    end type 

    type, public :: FlatTree 
        !! Flattened calculation tree for pointer-safe copying 

        type(CalculationKernel) ,allocatable ,dimension(:) :: kernels 
        type(IntArray) ,allocatable ,dimension(:) :: children 
        integer(ik) ,allocatable ,dimension(:) :: parent

    end type FlatTree

    type ,public ,extends(Object) :: CalculationNode
        !! Node class for the abstract calculation left child/right sibling tree. Each node has a reference to its leftmost child,
        !! its sibling to the immediate right, and to its parent. 
        !!
        !! Each node has a binary mode that is either additive or multiplicative, with 
        !! multiplicative being the default. This mode is applied to the results of the nodes children, either multiplying or adding them.
        !! 
        !! Each node also has a constant (defaulting to 1 in multiplicative and 0 in additive mode), which is applied to the node's result.
        !!
        !! An optional unary operation can be associated with each node, which is then applied to the node's result after the constant.
        !! The unary operation can be parameterized with real, integer, or logical parameters.
        !!
        !! Finally, if a node doesn't have a child (it is a leaf), it must have a nonzero index associated with a passed RealArray(:)
        !! list of variables, with the associated values then treated as the missing child's result.  

        type(CalculationNode) ,pointer :: leftChild => null()
        type(CalculationNode) ,pointer :: rightSibling => null()
        type(CalculationNode) ,pointer :: parent => null()

        type(CalculationKernel) :: kernel 

        procedure(realArrayFunctionGenParam) ,pointer ,nopass :: unaryTransform => null() !! Optional unary transformation

        contains

        procedure ,public :: addChild 
        procedure ,public :: init => initNode
        procedure ,public :: evaluate => evaluateNode
        procedure ,public :: destroy => destroyNode

    end type CalculationNode
!-----------------------------------------------------------------------------------------------------------------------------------
    type ,extends(Object) ,public:: CalculationTree

        type(CalculationNode), pointer :: root => null()

        contains

        procedure ,public :: init => initTree
        procedure ,public :: evaluate => evaluateTree
        procedure ,public :: flatten => flattenTree
        procedure ,public :: initFromFlatTree 

        final :: finalizeCalculationTree

    end type CalculationTree
!-----------------------------------------------------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------------------------------------------------
    interface 
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

        end subroutine addChild
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module recursive function evaluateNode(this,inputArray) result(res)
            !! Recursively evaluate nodes, using the inputArray variables for leaf values

            class(CalculationNode)        ,intent(in) :: this
            type(RealArray) ,dimension(:) ,intent(in) :: inputArray
            real(rk) ,allocatable ,dimension(:)       :: res
 
        end function evaluateNode
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module function evaluateTree(this,inputArray) result(res)
            !! Call tree's root node evaluate

            class(CalculationTree)        ,intent(in) :: this
            type(RealArray) ,dimension(:) ,intent(in) :: inputArray
            real(rk) ,allocatable ,dimension(:)       :: res
 
        end function evaluateTree
!-----------------------------------------------------------------------------------------------------------------------------------
        module function flattenTree(this) result(res)
            !! Flatten tree into FlatTree object

            class(CalculationTree)        ,intent(in) :: this
            type(FlatTree)                            :: res
 
        end function flattenTree
!-----------------------------------------------------------------------------------------------------------------------------------
        module subroutine initFromFlatTree(this,fTree)
            !! Calculation tree initialization routine using a FlatTree object

            class(CalculationTree)           ,intent(inout)  :: this
            type(FlatTree)                   ,intent(in)     :: fTree

        end subroutine initFromFlatTree
!-----------------------------------------------------------------------------------------------------------------------------------
        pure module recursive subroutine destroyNode(this)

            class(CalculationNode),intent(inout) :: this

        end subroutine destroyNode
!-----------------------------------------------------------------------------------------------------------------------------------
        elemental module subroutine finalizeCalculationTree(this) 

            type(CalculationTree) ,intent(inout) :: this

        end subroutine finalizeCalculationTree 
!-----------------------------------------------------------------------------------------------------------------------------------
    end interface
!-----------------------------------------------------------------------------------------------------------------------------------
 end module calculation_tree_class
!-----------------------------------------------------------------------------------------------------------------------------------
 