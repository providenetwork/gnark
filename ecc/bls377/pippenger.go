package bls377

// git branch experimental/yet-another-pippenger
// Gus's experiments on Pippenger's algorithm

import (
	"github.com/consensys/gnark/ecc/bls377/fr"
	"github.com/consensys/gnark/internal/debug"
)

// Gus549 implements the "Pippenger approach" from Section 4 of
// https://eprint.iacr.org/2012/549.pdf
func (p *G1Jac) Gus549(curve *Curve, points []G1Affine, scalars []fr.Element, c int) *G1Jac {
	// const c int = 4                        // scalars partitioned into c-bit radixes, must divide 64
	t := fr.ElementLimbs * 64 / c        // number of c-bit radixes in a scalar
	selectorMask := uint64((1 << c) - 1) // low c bits are 1

	debug.Assert(64%c == 0) // see TODO below
	debug.Assert(c < 31)

	buckets := make([]G1Jac, 1<<(c-1)) // only 2**(c-1) buckets instead of 2**c - 1

	// preprocess scalars into NAF form
	// TODO this is probably slower than it needs to be
	scalarsNAF := make([][]int32, len(points))
	for i := 0; i < len(points); i++ {
		scalarsNAF[i] = make([]int32, t)
		for j := 0; j < t; j++ {
			jc := j * c
			selectorIndex := jc / 64
			selectorShift := jc - (selectorIndex * 64)
			scalarsNAF[i][j] += int32((scalars[i][selectorIndex] & (selectorMask << selectorShift)) >> selectorShift)

			// ensure all scalars are at most 2**(c-1) in absolute value
			if scalarsNAF[i][j] > (1 << (c - 1)) {
				scalarsNAF[i][j] -= 1 << c
				scalarsNAF[i][j+1]++ // TODO assume this never happens for j==t-1
			}
		}
	}

	// notation: i ranges over points, scalars
	// notation: j ranges over c-bit radixes in a scalar
	// notation: s[i][j] := the jth c-bit radix of scalar[i]
	//
	// for each j:
	//   compute total[j] := 2^(j*c) * ( sum over i: s[i][j] * points[i] )
	// result p := ( sum over j: total[j] )

	for j := t - 1; j >= 0; j-- {

		// initialize 2^c - 1 buckets
		for k := 0; k < len(buckets); k++ {
			buckets[k].Set(&curve.g1Infinity)
		}

		// place points into buckets based on their selector
		for i := 0; i < len(points); i++ {
			if scalarsNAF[i][j] > 0 {
				buckets[scalarsNAF[i][j]-1].AddMixed(&points[i])
			} else if scalarsNAF[i][j] < 0 {
				var negPoint G1Affine
				negPoint.Neg(&points[i])
				buckets[-scalarsNAF[i][j]-1].AddMixed(&negPoint)
			}
		}

		// accumulate buckets into totalj
		var sumj, totalj G1Jac
		sumj.Set(&curve.g1Infinity)
		totalj.Set(&curve.g1Infinity)
		for k := len(buckets) - 1; k >= 0; k-- {
			sumj.Add(curve, &buckets[k])
			totalj.Add(curve, &sumj)
		}

		// accumulate totalj into result
		// if this is not the first iteration
		// then double p c times first
		if j == t-1 {
			p.Set(&totalj)
		} else {
			for l := 0; l < c; l++ {
				p.Double()
			}
			p.Add(curve, &totalj)
		}
	}

	return p
}

// Gus549 implements the "Pippenger approach" from Section 4 of
// https://eprint.iacr.org/2012/549.pdf
func (p *G1Jac) Gus549Old(curve *Curve, points []G1Affine, scalars []fr.Element, c int) *G1Jac {
	// const c int = 4                        // scalars partitioned into c-bit radixes, must divide 64
	t := fr.ElementLimbs * 64 / c        // number of c-bit radixes in a scalar
	selectorMask := uint64((1 << c) - 1) // low c bits are 1

	debug.Assert(64%c == 0) // see TODO below

	buckets := make([]G1Jac, (1<<c)-1)

	// notation: i ranges over points, scalars
	// notation: j ranges over c-bit radixes in a scalar
	// notation: s[i][j] := the jth c-bit radix of scalar[i]
	//
	// for each j:
	//   compute total[j] := 2^(j*c) * ( sum over i: s[i][j] * points[i] )
	// result p := ( sum over j: total[j] )

	for j := t - 1; j >= 0; j-- {

		// initialize 2^c - 1 buckets
		for k := 0; k < len(buckets); k++ {
			buckets[k].Set(&curve.g1Infinity)
		}

		// place points into buckets based on their selector
		jc := j * c
		selectorIndex := jc / 64
		selectorShift := jc - (selectorIndex * 64)
		for i := 0; i < len(points); i++ {

			// TODO: if c does not divide 64
			// then a c-bit radix might straddle two limbs of a scalar
			// -> need to fix this code
			selector := (scalars[i][selectorIndex] & (selectorMask << selectorShift)) >> selectorShift

			if selector != 0 {
				// TODO: don't add the first point.  Instead, check if the bucket is the zero point, and if so use Set instead of Add
				buckets[selector-1].AddMixed(&points[i])
			}
		}

		// accumulate buckets into totalj
		var sumj, totalj G1Jac
		sumj.Set(&curve.g1Infinity)
		totalj.Set(&curve.g1Infinity)
		for k := len(buckets) - 1; k >= 0; k-- {
			sumj.Add(curve, &buckets[k])
			totalj.Add(curve, &sumj)
		}

		// accumulate totalj into result
		// if this is not the first iteration
		// then double p c times first
		if j == t-1 {
			p.Set(&totalj)
		} else {
			for l := 0; l < c; l++ {
				p.Double()
			}
			p.Add(curve, &totalj)
		}
	}

	return p
}
