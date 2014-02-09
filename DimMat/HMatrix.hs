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
{- | This module exposes linear algebra backed by the hmatrix package.
-}
module DimMat.HMatrix (
  Matrix,
  Vector,
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
  Matrix :: (H.Container H.Matrix e, H.Field e, N.NumType (TListLength rs), N.NumType (TListLength cs))
         => H.Matrix e -> Matrix (MatrixShape g rs cs) e

-- define a type for vectors
-- Should be at kind VectorShape -> * -> *
data Vector shape e where
  Vector :: (H.Container H.Vector e, H.Field e)
         => H.Vector e -> Vector shape e


zeroes :: forall g rs cs e.(H.Field e, N.NumType (TListLength rs), N.NumType (TListLength cs)) => Matrix (MatrixShape g rs cs) e
zeroes = Matrix (H.konst 0
        (N.toNum (undefined :: TListLength rs),
         N.toNum (undefined :: TListLength cs)))
