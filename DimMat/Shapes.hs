{-# LANGUAGE CPP #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE OverlappingInstances #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE RankNTypes #-}
{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{- | This module only exposes the types and type functions neccessary to express linear algebra, it doesn't actually implement term-level linear algebra.
-}
module DimMat.Shapes (
  MatrixShape,
  Product,
  ShapeProduct,
  ShapeTranspose,
  ShapeInverse,
  ShapeDeterminant,
  ShapeDimensionless,
  ShapeRows,
  ShapeCols,
  ) where
    
import Data.Void (Void)
import GHC.Exts (Constraint)
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import qualified Numeric.NumType.TF as N

-- define a data kind for matrix shapes
-- a matrix shape is a single global dimension, an n-1 list of row dimesnions, and an m-1 list of column dimensions
-- Should be at kind Dim -> [Dim] -> [Dim] -> Dim
data MatrixShape :: * -> [*] -> [*] -> *


-- Define the type of matrix products.
-- Both arguments should be kinded as MatrixShape, as should result kind of ShapeProduct
class Product ldims rdims where
  type ShapeProduct ldims rdims :: *

instance (lcols ~ MapInv rrows) => Product (MatrixShape g1 rs1 cs1) (MatrixShape g2 rs2 cs2) where
  type ShapeProduct (MatrixShape g1 rs1 cs1) (MatrixShape g2 rs2 cs2) = MatrixShape (Mul g1 g2) rs1 cs2


-- Define the shape of matrix transposition.
-- Should be at kind MatrixShape -> MatrixShape
-- Should be a closed type family.
type family ShapeTranspose (mdims :: *) :: *

-- The shape of matrix transposition for shapes expressed in this form is simple, just flip the rows and columns.
type instance ShapeTranspose (MatrixShape g rs cs) = MatrixShape g cs rs


-- Define the shape of matrix inversion.
-- Should be at kind MatrixShape -> MatrixShape
-- Should be a closed type family.
type family ShapeInverse (mdims :: *) :: *

type instance ShapeInverse (MatrixShape g rs cs) = MatrixShape (Inv g) (MapInv cs) (MapInv rs)


-- Define the type of a matrix determinant.
-- (This is defined even where the determinant doesn't exist (non-square shapes), but no problem arises because you can't call the actual term-level determinant method at non-square shapes.
-- Should be at kind MatrixShape -> Dim
-- Should be a closed type family.
type family ShapeDeterminant (shape :: *) :: *

type instance ShapeDeterminant (MatrixShape g rs cs) = Mul g (Mul (DimProduct rs) (DimProduct cs))


-- Define the type of a matrix the same size, but with dimensions stripped.
-- Should be at kind MatrixShape -> MatrixShape
-- Should be a closed type family.
type family ShapeDimensionless (shape :: *) :: *

type instance ShapeDimensionless (MatrixShape g rs cs) = MatrixShape DOne (MapConstOne cs) (MapConstOne cs)


-- Define the type-level number of rows in a matrix.
-- Should be at kind MatrixShape -> Nat
-- Should be a closed type family.
type family ShapeRows (shape :: *) :: *

type instance ShapeRows (MatrixShape g rs cs) = N.S (TListLength rs)


-- Define the type-level number of columns in a matrix.
-- Should be at kind MatrixShape -> Nat
-- Should be a closed type family.
type family ShapeCols (shape :: *) :: *

type instance ShapeCols (MatrixShape g rs cs) = N.S (TListLength cs)






-- Shorthand for the inverse of a dimension.
-- Should be at kind Dim -> Dim
type Inv d = Div DOne d


-- Invert all dimensions in a list of dimensions.
-- Should be at kind [Dim] -> [Dim]
-- Should be a closed type family.
type family MapInv (dims :: [*]) :: [*]

type instance MapInv '[] = '[]
type instance MapInv (x ': xs) = (Inv x) ': (MapInv xs)


-- Convert all dimensions in a list of dimensions to dimensionless.
-- Should be at kind [Dim] -> [Dim]
-- Should be a closed type family.s
type family MapConstOne (dims :: [*]) :: [*]

type instance MapConstOne '[] = '[]
type instance MapConstOne (x ': xs) = DOne ': (MapConstOne xs)


-- Define the product of a list of dimensions.
-- Should be at kind [Dim] -> Dim
-- Should be a closed type family.
type family DimProduct (dims :: [*]) :: *

type instance DimProduct '[] = DOne
type instance DimProduct (a ': as) = Mul a (DimProduct as)


-- Get the length of a type-level list.
-- Should be at kind [k] -> Nat
-- Should be a closed type family.
type family TListLength (xs :: [*]) :: *

type instance TListLength '[] = N.Z
type instance TListLength (x ': xs) = N.S (TListLength xs)


-- Get a specified, zero-indexed element from a type-level list.
-- Should be at kind [k] -> k (except what about the Void case?)
-- Should be a closed type family.
type family TElementAt (xs :: [*]) (n :: *) :: *

type instance TElementAt '[] n = Void
type instance TElementAt (a ': as) N.Z = a
type instance TElementAt (a ': as) (N.S n) = TElementAt as n
