{-# LANGUAGE ScopedTypeVariables #-}

{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE PolyKinds #-}
{-# LANGUAGE ConstraintKinds #-}

{-# LANGUAGE MultiParamTypeClasses #-}

{-# LANGUAGE TypeFamilies #-}
{-# LANGUAGE TypeOperators #-}

{-# LANGUAGE UndecidableInstances #-}

{-# LANGUAGE GADTs #-}
{-# LANGUAGE FlexibleContexts #-}

{-# LANGUAGE StandaloneDeriving #-}
{- | This module exposes linear algebra backed by the hmatrix package.
-}
module DimMat.HMatrix (
  Matrix,
  Vector,
  zeroes,
  zeroesV,
  add,
  sub,
  rank,
  rows,
  cols,
  ) where

import Data.Void (Void)
import GHC.Exts (Constraint)
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import qualified Numeric.NumType.TF as N
import Data.List.NonEmpty (NonEmpty(..))

import qualified Numeric.LinearAlgebra as H
import qualified Numeric.LinearAlgebra.LAPACK as H

import DimMat.Shapes

-- define a type for matrices
-- Should be at kind MatrixShape -> * -> *
data Matrix shape e where
  Matrix :: (H.Container H.Matrix e, H.Field e, N.NumType (ShapeRows shape), N.NumType (ShapeCols shape))
         => H.Matrix e -> Matrix shape e

-- very basic show instance, doesn't include dimensions at all
deriving instance (Show e) => Show (Matrix shape e)

-- define a type for vectors
-- Should be at kind VectorShape -> * -> *
data Vector shape e where
  Vector :: (H.Container H.Vector e, H.Field e)
         => H.Vector e -> Vector shape e

-- very basic show instance, doesn't include dimensions at all
deriving instance (Show e) => Show (Vector shape e)


-- a matrix of zeroes
zeroes :: forall shape e.(H.Field e, N.NumType (ShapeRows shape), N.NumType (ShapeCols shape)) => Matrix shape e
zeroes = Matrix (H.konst 2
        (N.toNum (undefined :: ShapeRows shape),
         N.toNum (undefined :: ShapeCols shape)))

-- a vector of zeroes
zeroesV :: forall shape e.(H.Field e, N.NumType (VectorLength shape)) => Vector shape e
zeroesV = Vector (H.constant 0 (N.toNum (undefined :: VectorLength shape)))

-- some simple matrix functions
-- could implement with liftMatrix, but that loses the GADT evidence for H.Field e
rank :: Matrix shape e -> Int
rank (Matrix m) = H.rank m

rows :: Matrix shape e -> Int
rows = liftMatrix H.rows

cols :: Matrix shape e -> Int
cols = liftMatrix H.cols

-- addition and subtraction
class DimensionalContainer c where
  add :: forall shape e.(H.Field e) => c shape e -> c shape e -> c shape e
  sub :: forall shape e.(H.Field e) => c shape e -> c shape e -> c shape e

instance DimensionalContainer Matrix where
  add = liftMatrix2 (H.add)
  sub = liftMatrix2 (H.sub)

instance DimensionalContainer Vector where
  add = liftVector2 (H.add)
  sub = liftVector2 (H.sub)

liftMatrix :: (m ~ Matrix shape e,
               h ~ H.Matrix e) => (h -> t) -> m -> t
liftMatrix f (Matrix a) = f a

liftMatrix1 :: (m ~ Matrix shape e,
                h ~ H.Matrix e) => (h -> h) -> m -> m
liftMatrix1 f (Matrix a) = Matrix (f a)

liftMatrix2 :: (m ~ Matrix shape e,
                h ~ H.Matrix e) => (h -> h -> h) -> m -> m -> m
liftMatrix2 f (Matrix a) (Matrix b) = Matrix (f a b)

liftVector :: (v ~ Vector shape e,
               h ~ H.Vector e) => (h -> t) -> v -> t
liftVector f (Vector a) = f a

liftVector1 :: (v ~ Vector shape e,
                h ~ H.Vector e) => (h -> h) -> v -> v
liftVector1 f (Vector a) = Vector (f a)

liftVector2 :: (v ~ Vector shape e,
                h ~ H.Vector e) => (h -> h -> h) -> v -> v -> v
liftVector2 f (Vector a) (Vector b) = Vector (f a b)
