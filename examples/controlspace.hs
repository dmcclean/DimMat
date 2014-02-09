{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE MultiParamTypeClasses #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE NoMonomorphismRestriction #-}
{-# LANGUAGE ConstraintKinds #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE UndecidableInstances #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE StandaloneDeriving #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE QuasiQuotes #-}
{-# LANGUAGE GADTs #-}
{-# LANGUAGE DataKinds #-}
{-# LANGUAGE KindSignatures #-}
{-# LANGUAGE TypeFamilies #-}
module T1 where
import GHC.Exts (Constraint)
import DimMat
import Numeric.Units.Dimensional.TF.Prelude
import Numeric.Units.Dimensional.TF
import qualified Prelude as P
import Text.PrettyPrint.ANSI.Leijen hiding ((<>))


import qualified Numeric.LinearAlgebra as H

{- Example from http://ctms.engin.umich.edu/CTMS/index.php?example=InvertedPendulum&section=ControlStateSpace -}

-- not really stated on that page, but if you go back a couple of pages in their derivation
-- you can see that the type of u is a 1x1 matrix whose sole element is a force

massOfCart = (0.5 :: Double) *~ (kilo gram)
massOfPendulum = (0.2 :: Double) *~ (kilo gram)
coefficientOfFrictionForCart = (0.1 :: Double) *~ (newton / (meter / second))
lengthToPendulumCenterOfMass = (0.3 :: Double) *~ meter
massMomentOfInertiaOfPendulum = (0.006 :: Double) *~ (kilo gram * meter^pos2)
g = (9.8 :: Double) *~ (meter / second^pos2)

p = massMomentOfInertiaOfPendulum*(massOfCart+massOfPendulum)+(massOfCart*massOfPendulum*lengthToPendulumCenterOfMass*lengthToPendulumCenterOfMass)

a22 = negate (massMomentOfInertiaOfPendulum+massOfPendulum * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass) * coefficientOfFrictionForCart / p
a23 = (massOfPendulum * massOfPendulum * g * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass) / p
a42 = negate (massOfPendulum * lengthToPendulumCenterOfMass * coefficientOfFrictionForCart) / p
a43 = massOfPendulum * g * lengthToPendulumCenterOfMass*(massOfCart + massOfPendulum)/p

b21 = (massMomentOfInertiaOfPendulum + (massOfPendulum * lengthToPendulumCenterOfMass * lengthToPendulumCenterOfMass)) / p
b41 = massOfPendulum * lengthToPendulumCenterOfMass / p

-- example state value
x = [matD| 1.0 *~ meter;
           0.2 *~ (meter / second);
           _0 :: Dimensionless Double;
           0.1 *~ (second^neg1) |]

dx = scale (_1 / (1 *~ second)) x

-- example control input
u = [matD| (0 :: Double) *~ newton |]

--xDot = (a `multiply` x) `add` (b `multiply` u)

y = [matD| 1 *~ meter; _0 :: Dimensionless Double  |]

poles sys = let ContinuousLiSystem { a'' = a } = sys
             in eigenvalues a

-- this ludicrous context doesn't actually mean anything at all, and is always met for well-kinded systems
discretizeZeroOrderHold :: (H.Field e,
                            SameLength (MapDiv (Head xs) xs)
                                       (MapMul (Head xs) (MapRecip xs)),
                            SameLength (MapMul (Head xs) (MapRecip xs))
                                       (MapDiv (Head xs) xs),
                            AreRecipsList (MapDiv (Head xs) xs)
                                          (MapMul (Head xs) (MapRecip xs)),
                            MapMultEq iv (MapDiv (Head xs) (MapDiv iv xs))
                                         (MapDiv (Head xs) xs),
                            (MapDiv (Head us) xs) ~ (t1 ': t2),
                            HNat2Integral (HLength t2),
                            (MapMul (Head us) (MapRecip us)) ~ (DOne ': t3),
                            HNat2Integral (HLength t3)
                           ) => ContinuousLiSystem iv xs ys us e -> Quantity iv e -> DiscreteLiSystem iv xs ys us e
discretizeZeroOrderHold sys t = let ac = a'' sys
                                    bc = b'' sys
                                    ad = expm (scale t (a'' sys))
                                    --bd = (pinv ac) <> (ad `add` (scale (negate _1) ident)) <> bc
                                    bd = zeroes
                                 in DiscreteLiSystem { t''' = t, a''' = ad, b''' = bd, c''' = c'' sys, d''' = d'' sys }

ctrb = let ContinuousLiSystem { a'' = a, b'' = b, c'' = c, d'' = d } = pendulum
        in [blockD| b, a <> b, a <> a <> b, a <> a <> a <> b |]

outputctrb = let ContinuousLiSystem { a'' = a, b'' = b, c'' = c, d'' = d } = pendulum
                 cb = c <> b
                 cab = c <> a <> b
                 caab = c <> a <> a <> b
                 caaab = c <> a <> a <> b
              in [blockD| cb, cab, caab, caaab, d |]

obsv = let ContinuousLiSystem { a'' = a, b'' = b, c'' = c, d'' = d } = pendulum
        in [blockD| c;
                    c <> a;
                    c <> a <> a;
                    c <> a <> a <> a |]

type family DivideVectors (to :: [*]) (from :: [*]) :: [[*]]
type instance DivideVectors ts fs = [ MapDiv (Head fs) ts, MapMul (Head fs) (MapRecip fs) ]

data ContinuousLiSystem (iv :: *) (xs :: [*]) (ys :: [*]) (us :: [*]) e = ContinuousLiSystem
                                                                          {
                                                                            a'' :: DimMat (DivideVectors (MapDiv iv xs) xs) e,
                                                                            b'' :: DimMat (DivideVectors (MapDiv iv xs) us) e,
                                                                            c'' :: DimMat (DivideVectors ys xs) e,
                                                                            d'' :: DimMat (DivideVectors ys us) e
                                                                          }

{-
deriving instance (Show e, 
                   PPUnits (DivideVectors (MapDiv iv xs) xs),
                   PPUnits (DivideVectors (MapDiv iv xs) us),
                   PPUnits (DivideVectors ys xs),
                   PPUnits (DivideVectors ys us)) => Show (ContinuousLiSystem iv xs ys us e)
-}

type ContinuousLtiSystem = ContinuousLiSystem DTime

data DiscreteLiSystem (iv :: *) (xs :: [*]) (ys :: [*]) (us :: [*]) e = DiscreteLiSystem
                                                                        {
                                                                          t''' :: Quantity iv e,
                                                                          a''' :: DimMat (DivideVectors xs xs) e,
                                                                          b''' :: DimMat (DivideVectors xs us) e,
                                                                          c''' :: DimMat (DivideVectors ys xs) e,
                                                                          d''' :: DimMat (DivideVectors ys us) e
                                                                        }

{-
deriving instance (Show e,
                   Show iv,
                   PPUnits (DivideVectors xs xs),
                   PPUnits (DivideVectors xs us),
                   PPUnits (DivideVectors ys xs),
                   PPUnits (DivideVectors ys us)) => Show (DiscreteLiSystem iv xs ys us e)
-}

type DiscreteLtiSystem = DiscreteLiSystem DTime

type SimpleExampleSystem = ContinuousLiSystem DOne
    '[DOne]
    '[DOne]
    '[DOne]
    Double

type ExampleSystem = ContinuousLtiSystem
    '[DLength, DVelocity, DPlaneAngle, DAngularVelocity]
    '[DLength, DPlaneAngle]
    '[DForce]
    Double

simple = ContinuousLiSystem {
           a'' = zeroes,
           b'' = zeroes,
           c'' = ident,
           d'' = zeroes
         } :: SimpleExampleSystem

pendulum = ContinuousLiSystem {
           a'' = [matD| _0, _1, _0, _0;
                        _0, a22, a23, _0;
                        _0, _0, _0, _1;
                        _0, a42, a43, _0 |],
           b'' = [matD| _0;
                        b21;
                        _0;
                        b41 |],
           c'' = [matD| _1, _0, _0, _0;
                        _0, _0, _1, _0 |],
           d'' = zeroes
         } :: ExampleSystem

--evaluate :: LiSystem iv xs ys us e -> DimMat [xs,'[DOne]] e -> DimMat [us, '[DOne]] e -> (DimMat [(MapDiv iv xs), '[DOne]] e, DimMat [ys, '[DOne]] e)
evaluate sys x u = let
                      a = a'' sys
                      b = b'' sys
                      c = c'' sys
                      d = d'' sys
                      xDot = (a <> x) `add` (b <> u)
                      y = (c <> x) `add` (d <> u)
                    in (xDot, y)

evaluateD sys x u = let
                       a = a''' sys
                       b = b''' sys
                       c = c''' sys
                       d = d''' sys
                       x' = (a <> x) `add` (b <> u)
                       y = (c <> x) `add` (d <> u)
                     in (x', y)