# HiGarrote 1.0.1

## Bug fixes

* Remove the original make_PD function in lib.cpp.
* Use Matrix::nearPD instead.

# HiGarrote 1.0.2

## Bug fixes

* Remove the original is_PD function in lib.cpp.
* Use matrixcalc::is.positive.definite instead.

# HiGarrote 1.1.0

## Changes in version 1.0.3 (submitted 2025-04-20)

* Remove the original is_PD function in function: nnGarrote and use matrixcalc::is.positive.definite instead.
* Resolve the scaling issues happening in functions: HiGarrote and nnGarrote in HiGarrote.R. The new version ensures that variable selection is invariant to the response scaling.