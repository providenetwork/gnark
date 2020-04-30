package bls377

// git branch experimental/yet-another-pippenger
// Gus's experiments on Pippenger's algorithm

import (
	"fmt"
	"runtime"
	"testing"

	"github.com/consensys/gnark/ecc/bls377/fr"
)

func TestGus549(t *testing.T) {
	const fewPoints = 30
	const manyPoints = 3000
	cs := [...]int{4, 8, 16}

	curve := BLS377()

	var G G1Jac

	// mixer ensures that all the words of a fpElement are set
	var mixer fr.Element
	mixer.SetString("7716837800905789770901243404444209691916730933998574719964609384059111546487")

	samplePoints := make([]G1Affine, manyPoints)
	sampleScalars := make([]fr.Element, manyPoints)

	G.Set(&curve.g1Gen)

	for i := 1; i <= manyPoints; i++ {
		sampleScalars[i-1].SetUint64(uint64(i)).
			MulAssign(&mixer).
			FromMont()
		G.ToAffineFromJac(&samplePoints[i-1])
		G.Add(curve, &curve.g1Gen)
	}

	var finalBigScalar fr.Element
	var finalLotOfPoint G1Jac
	finalBigScalar.SetString("9004500500").MulAssign(&mixer).FromMont()
	finalLotOfPoint.ScalarMul(curve, &curve.g1Gen, finalBigScalar)

	var finalScalar fr.Element
	var finalPoint G1Jac
	finalScalar.SetString("9455").MulAssign(&mixer).FromMont()
	finalPoint.ScalarMul(curve, &curve.g1Gen, finalScalar)

	var testLotOfPoint, testPoint G1Jac
	for _, c := range cs {
		testLotOfPoint.Gus549(curve, samplePoints, sampleScalars, c)
		testPoint.Gus549(curve, samplePoints[:fewPoints], sampleScalars[:fewPoints], c)

		if !finalLotOfPoint.Equal(&testLotOfPoint) {
			t.Fatal("error multi (>50 points) exp")
		}
		if !finalPoint.Equal(&testPoint) {
			t.Fatal("error multi <=50 points) exp")
		}

		testLotOfPoint.Gus549Old(curve, samplePoints, sampleScalars, c)
		testPoint.Gus549Old(curve, samplePoints[:fewPoints], sampleScalars[:fewPoints], c)

		if !finalLotOfPoint.Equal(&testLotOfPoint) {
			t.Fatal("error multi (>50 points) exp")
		}
		if !finalPoint.Equal(&testPoint) {
			t.Fatal("error multi <=50 points) exp")
		}

	}
}

func testPointsG1MultiExp(n int) (points []G1Jac, scalars []fr.Element) {

	curve := BLS377()

	// points
	points = make([]G1Jac, n)
	points[0].Set(&curve.g1Gen)
	points[1].Set(&points[0]).Double() // can't call p.Add(a) when p equals a
	for i := 2; i < len(points); i++ {
		points[i].Set(&points[i-1]).Add(curve, &points[0]) // points[i] = i*g1Gen
	}

	// scalars
	// non-Montgomery form
	// cardinality of G1 is the fr modulus, so scalars should be fr.Elements
	// non-Montgomery form
	scalars = make([]fr.Element, n)

	// To ensure a diverse selection of scalars that use all words of an fr.Element,
	// each scalar should be a power of a large generator of fr.
	var scalarGenMont fr.Element
	scalarGenMont.SetString("7716837800905789770901243404444209691916730933998574719964609384059111546487")
	scalars[0].Set(&scalarGenMont).FromMont()

	var curScalarMont fr.Element // Montgomery form
	curScalarMont.Set(&scalarGenMont)
	for i := 1; i < len(scalars); i++ {
		curScalarMont.MulAssign(&scalarGenMont) // scalars[i] = scalars[0]^i
		scalars[i].Set(&curScalarMont).FromMont()
	}

	return points, scalars
}

//--------------------//
//     benches		  //
//--------------------//

func BenchmarkGus549(b *testing.B) {

	fmt.Println("GOMAXPROCS was", runtime.GOMAXPROCS(1))

	curve := BLS377()
	numPoints := [...]int{4e4, 1e5, 1e6}
	var exp G1Jac
	cs := [...]int{8, 16}

	for j := range numPoints {
		points, scalars := testPointsG1MultiExp(numPoints[j])

		// MultiExp takes affine points, not Jacobian points
		pointsAffine := make([]G1Affine, len(points))
		for k := range pointsAffine {
			points[k].ToAffineFromJac(&pointsAffine[k])
		}

		for _, c := range cs {
			b.Run(fmt.Sprintf("%d-Gus549-%d", numPoints[j], c), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					exp.Gus549(curve, pointsAffine, scalars, c)
				}
			})
			b.Run(fmt.Sprintf("%d-Gus549Old-%d", numPoints[j], c), func(b *testing.B) {
				for i := 0; i < b.N; i++ {
					exp.Gus549Old(curve, pointsAffine, scalars, c)
				}
			})
		}

		b.Run(fmt.Sprintf("%d-MultiExp", numPoints[j]), func(b *testing.B) {
			for i := 0; i < b.N; i++ {
				<-exp.MultiExp(curve, pointsAffine, scalars)
			}
		})

		// b.Run(fmt.Sprintf("%d-multiExp", numPoints[j]), func(b *testing.B) {
		// 	for i := 0; i < b.N; i++ {
		// 		exp.multiExp(curve, points, scalars)
		// 	}
		// })
	}
}
